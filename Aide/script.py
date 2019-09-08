#! /usr/bin/env python3
import math
import sys

def calc_distance3D(A, B):
    """
    Fonction qui permet de calculer la distance 3D entre deux listes de coordonnées A et B

    :param A: A[0] = x, A[1] = y, A[2] = z
    :param B: B[0] = x, B[1] = y, B[2] = z
    :return: distance entre deux points 3D
    """
    d = math.sqrt((B[0]-A[0])**2+(B[1]-A[1])**2+(B[2]-A[2])**2)
    return d


def trouve_charge(resi, atome, file):
    """
    Fonction qui va aller chercher la charge d'un atome donné d'un résidus
    donné dans le fichier de configuration de type CHARMM donné
    :param resi: Nom du résidus voulus
    :param atome: Nom de l'atome voulus
    :param file: Nom du fichier de configuration
    :return: Charge de l'atome voulus
    """
    bon_res = 0
    with open(file, 'r') as fichier:
        for ligne in fichier:
            if ligne[5:8] == resi:
                bon_res = 1
            if bon_res == 1 and ligne[0:10] == "ATOM " + atome + " "*(5-len(atome)):
                return "".join(ligne[17:23])


def trouve_ASN_OD1_ND2_GLN_OE1_NE2(nom_pdb):
    """
    Fonction qui va Trouver tout les atomes OD1 et ND2 de toutes les ASN
    ainsi que tout les atomes OE1 et NE2 de tout les GLN
    :param nom_pdb: Nom du fichier pdb a analyser
    :return: 4 listes pour chaque atomes (un indice = 1 ligne du fichier pdb)
    """
    ligne_ASN_OD1 = []
    ligne_ASN_ND2 = []
    ligne_GLN_OE1 = []
    ligne_GLN_NE2 = []
    with open(nom_pdb, 'r') as fichier:
        ligne_fichier = fichier.readline()
        while ligne_fichier != "":
            ligne_fichier = fichier.readline()
            if "ASN" == ligne_fichier[17:20].strip() \
                    and "ATOM" == ligne_fichier[0:6].strip() \
                    and "OD1" == ligne_fichier[13:16].strip():
                ligne_ASN_OD1.append(ligne_fichier)
            if "ASN" == ligne_fichier[17:20].strip() \
                    and "ATOM" == ligne_fichier[0:6].strip() \
                    and "ND2" == ligne_fichier[13:16].strip():
                ligne_ASN_ND2.append(ligne_fichier)

            if "GLN" == ligne_fichier[17:20].strip() \
                    and "ATOM" == ligne_fichier[0:6].strip() \
                    and "OE1" == ligne_fichier[13:16].strip():
                ligne_GLN_OE1.append(ligne_fichier)
            if "GLN" == ligne_fichier[17:20].strip() \
                    and "ATOM" == ligne_fichier[0:6].strip() \
                    and "NE2" == ligne_fichier[13:16].strip():
                ligne_GLN_NE2.append(ligne_fichier)
        return ligne_ASN_OD1, ligne_ASN_ND2, ligne_GLN_OE1, ligne_GLN_NE2


def trouve_atom_N_O(nom_pdb):
    """
    Fonction qui trouve tous les atomes N ou O des protéines présent dans le fichier pdb
    :param nom_pdb: Nom du fichier pdb ou chercher
    :return: Liste de tout les atomes N ou O présent
            Liste[i][0] = coord x
            Liste[i][1] = coord y
            Liste[i][2] = coord z
            Liste[i][3] = Nom du résidus
            Liste[i][4] = Nom de l'atomes
    """
    liste_coord_6A = []
    with open(nom_pdb, 'r') as fichier:
        ligne_fichier = fichier.readline()
        while ligne_fichier != "":
            ligne_fichier = fichier.readline()
            if "ATOM" == ligne_fichier[0:6].strip()\
                and "N" == ligne_fichier[13] \
                and "PRO" == ligne_fichier[72:75].strip():
                coord_6A = []
                x = float(ligne_fichier[30:38].strip())
                coord_6A.append(x)
                y = float(ligne_fichier[38:46].strip())
                coord_6A.append(y)
                z = float(ligne_fichier[46:54].strip())
                coord_6A.append(z)
                residue = ligne_fichier[17:20].strip()
                coord_6A.append(residue)
                atom = ligne_fichier[13:16].strip()
                coord_6A.append(atom)

                liste_coord_6A.append(coord_6A)

            if "ATOM" == ligne_fichier[0:6].strip() \
                and "O" == ligne_fichier[13] \
                and "PRO" == ligne_fichier[72:75].strip():
                coord_6A = []
                x = float(ligne_fichier[30:38].strip())
                coord_6A.append(x)
                y = float(ligne_fichier[38:46].strip())
                coord_6A.append(y)
                z = float(ligne_fichier[46:54].strip())
                coord_6A.append(z)
                residue = ligne_fichier[17:20].strip()
                coord_6A.append(residue)
                atom = ligne_fichier[13:16].strip()
                coord_6A.append(atom)

                liste_coord_6A.append(coord_6A)

        return liste_coord_6A


def liste_to_coord(liste):
    """
    Fonction qui va extraire les coordonnées d'une ligne pdb
    Ici d'une liste de ligne pdb
    :param liste: liste de ligne pdb
    :return: une liste de coordonnées
    """
    liste_coord = []
    for i in range(len(liste)):
        coord = []
        x = float(liste[i][30:38].strip())
        coord.append(x)
        y = float(liste[i][38:46].strip())
        coord.append(y)
        z = float(liste[i][46:54].strip())
        coord.append(z)
        liste_coord.append(coord)
    return liste_coord


def fct_calcul_energie(nom_residue, nom_atome_O, liste_residue_O, nom_atome_N, liste_residue_N, fichier_topo):
    """
    Fonction qui a pour but de ressortir les coordonnées des atomes dans l'ordre, c'est à dire si il y a FLIP ou pas.
    Etape 1 : Elle va calculer les distances entre chaque atomes et l'atome d'intérêt
    Etape 2 : Si la distance est <= a 6 angström, alors on calcul l'énergie électrostatique
    Etape 3 : Refais les étapes 1 et 2 pour le Flip
    Etape 4 : Evalue si FLIP ou non
    Etape 5 : Crée une liste de liste ou les coordonnées sont dans l'ordre

    :param nom_residue: Nom du résidue d'intérêt
    :param nom_atome_O: Nom de l'atome O du résidue d'intérêt
    :param liste_residue_O: Liste de l'atome O du résidue d'intérêt
    :param nom_atome_N: Nom de l'atome N du résidue d'intérêt
    :param liste_residue_N: Liste de l'atome N du résidue d'intérêt
    :param fichier_topo: Fichier de topologie
    :return: Liste de coordonnées dans l'ordre (FLIP ou non)
    """
    coord_dans_ordre = []
    flip = 0
    for i in range(len(liste_residue_O)):
        liste_coord_O = liste_to_coord(liste_residue_O)
        liste_coord_N = liste_to_coord(liste_residue_N)
        v_elec_O = 0
        v_elec_N = 0
        v_elec_O_flip = 0
        v_elec_N_flip = 0

        for j in range(len(liste_coord_N_O)):
            dist = calc_distance3D(liste_coord_O[i], liste_coord_N_O[j])
            if dist <= 6 and dist != 0:
                q_ASN_O = trouve_charge(nom_residue, nom_atome_O, fichier_topo)
                q_ASN_O = float(q_ASN_O)
                q_autre_O = trouve_charge(liste_coord_N_O[j][3], liste_coord_N_O[j][4], fichier_topo)
                q_autre_O = float(q_autre_O)
                v_elec = (q_ASN_O * q_autre_O) / (4 * math.pi * 1 * dist)
                v_elec_O += v_elec
                # print(v_elec_O)

            dist = calc_distance3D(liste_coord_N[i], liste_coord_N_O[j])
            if dist <= 6 and dist != 0:
                q_ASN_N = trouve_charge(nom_residue, nom_atome_N, fichier_topo)
                q_ASN_N = float(q_ASN_N)
                q_autre_N = trouve_charge(liste_coord_N_O[j][3], liste_coord_N_O[j][4], fichier_topo)
                q_autre_N = float(q_autre_N)
                v_elec = (q_ASN_N * q_autre_N) / (4 * math.pi * 1 * dist)
                v_elec_N += v_elec
                # print(v_elec_N)

            dist = calc_distance3D(liste_coord_O[i], liste_coord_N_O[j])
            if dist <= 6 and dist != 0:
                q_ASN_N = trouve_charge(nom_residue, nom_atome_N, fichier_topo)
                q_ASN_N = float(q_ASN_N)
                q_autre_N = trouve_charge(liste_coord_N_O[j][3], liste_coord_N_O[j][4], fichier_topo)
                q_autre_N = float(q_autre_N)
                v_elec = (q_ASN_N * q_autre_N) / (4 * math.pi * 1 * dist)
                v_elec_N_flip += v_elec
                # print(v_elec_N_flip)

            dist = calc_distance3D(liste_coord_N[i], liste_coord_N_O[j])
            if dist <= 6 and dist != 0:
                q_ASN_O = trouve_charge(nom_residue, nom_atome_O, fichier_topo)
                q_ASN_O = float(q_ASN_O)
                q_autre_O = trouve_charge(liste_coord_N_O[j][3], liste_coord_N_O[j][4], fichier_topo)
                q_autre_O = float(q_autre_O)
                v_elec = (q_ASN_O * q_autre_O) / (4 * math.pi * 1 * dist)
                v_elec_O_flip += v_elec
                # print(v_elec_O_flip)

        v_elec_sans_flip = v_elec_O + v_elec_N
        v_elec_avec_flip = v_elec_O_flip + v_elec_N_flip

        if v_elec_sans_flip < v_elec_avec_flip:
            coord_dans_ordre.append(liste_coord_O[i][0])
            coord_dans_ordre.append(liste_coord_O[i][1])
            coord_dans_ordre.append(liste_coord_O[i][2])
            coord_dans_ordre.append(liste_coord_N[i][0])
            coord_dans_ordre.append(liste_coord_N[i][1])
            coord_dans_ordre.append(liste_coord_N[i][2])
        elif v_elec_sans_flip > v_elec_avec_flip:
            flip += 1
            coord_dans_ordre.append(liste_coord_N[i][0])
            coord_dans_ordre.append(liste_coord_N[i][1])
            coord_dans_ordre.append(liste_coord_N[i][2])
            coord_dans_ordre.append(liste_coord_O[i][0])
            coord_dans_ordre.append(liste_coord_O[i][1])
            coord_dans_ordre.append(liste_coord_O[i][2])
        else:
            print("égalité")

    print("Il ya a eu {0} flip sur les {1} résidues {2}".format(flip, len(liste_residue_O), nom_residue))
    return coord_dans_ordre


def ecrire_fichier_pdb_avec_flip(fichier_pdb, fichier_ecriture, coord_dans_ordre_ASN, coord_dans_ordre_GLN):
    """
    Fonction qui va copier le fichier pdb d'entrée en faisant le flip des coordonnées si nécessaire
    :param fichier_pdb: Fichier pdb d'entrée
    :param fichier_ecriture: Fichier pdb de sortie
    :param coord_dans_ordre_ASN: Liste des coordonnées d'ASN dans l'ordre
    :param coord_dans_ordre_GLN: Liste des coordonnées de GLN dans l'ordre
    :return: Rien, écris juste le nouveu fichier pdb
    """
    nbr_ASN = 0
    nbr_GLN = 0

    with open(fichier_pdb, 'r') as fichier:
        with open(fichier_ecriture, 'w') as ecrit:
            for ligne in fichier:
                if ligne[17:20].strip() == "ASN"\
                        and "ATOM" == ligne[0:6].strip() \
                        and "OD1" == ligne[13:16].strip():
                    # or "ND2" == ligne[13:16].strip() ):
                    ecrit.write(ligne[0:30] + "{:8.3f}{:8.3f}{:8.3f}".format(coord_dans_ordre_ASN[nbr_ASN], coord_dans_ordre_ASN[nbr_ASN+1], coord_dans_ordre_ASN[nbr_ASN+2]) + ligne[54:80])
                    nbr_ASN += 3
                elif ligne[17:20].strip() == "ASN" \
                        and "ATOM" == ligne[0:6].strip() \
                        and "ND2" == ligne[13:16].strip():
                    ecrit.write(ligne[0:30] + "{:8.3f}{:8.3f}{:8.3f}".format(coord_dans_ordre_ASN[nbr_ASN], coord_dans_ordre_ASN[nbr_ASN+1], coord_dans_ordre_ASN[nbr_ASN+2]) + ligne[54:80])
                    nbr_ASN += 3
                elif ligne[17:20].strip() == "GLN"\
                        and "ATOM" == ligne[0:6].strip() \
                        and "OE1" == ligne[13:16].strip() :
                    # or :"NE2" == ligne[13:16].strip() ):
                    ecrit.write(ligne[0:30] + "{:8.3f}{:8.3f}{:8.3f}".format(coord_dans_ordre_GLN[nbr_GLN], coord_dans_ordre_GLN[nbr_GLN+1], coord_dans_ordre_GLN[nbr_GLN+2]) + ligne[54:80])
                    nbr_GLN += 3
                elif ligne[17:20].strip() == "GLN" \
                        and "ATOM" == ligne[0:6].strip() \
                        and "NE2" == ligne[13:16].strip():
                    ecrit.write(ligne[0:30] + "{:8.3f}{:8.3f}{:8.3f}".format(coord_dans_ordre_GLN[nbr_GLN], coord_dans_ordre_GLN[nbr_GLN+1], coord_dans_ordre_GLN[nbr_GLN+2]) + ligne[54:80])
                    nbr_GLN += 3
                else:
                    ecrit.write(ligne)


if len(sys.argv) < 3:
    # print("Il vous manque un argument, la bonne ligne de commande est \n python3 script.py file.pdb top_all36_prot.rtf")
    sys.exit("Il vous manque un argument, la bonne ligne de commande est \n \
     python3 script.py file.pdb top_all36_prot.rtf")



fichier_pdb = sys.argv[1]
fichier_topo = sys.argv[2]
fichier_ecriture = "test_ecriture_fichier.pdb"

ligne_ASN_OD1, ligne_ASN_ND2, ligne_GLN_OE1, ligne_GLN_NE2 = \
    trouve_ASN_OD1_ND2_GLN_OE1_NE2(fichier_pdb)

liste_coord_N_O = trouve_atom_N_O(fichier_pdb)

coord_dans_ordre_ASN = fct_calcul_energie("ASN", "OD1", ligne_ASN_OD1, "ND2", ligne_ASN_ND2, fichier_topo)
coord_dans_ordre_GLN = fct_calcul_energie("GLN", "OE1", ligne_GLN_OE1, "NE2", ligne_GLN_NE2, fichier_topo)

ecrire_fichier_pdb_avec_flip (fichier_pdb, fichier_ecriture, coord_dans_ordre_ASN, coord_dans_ordre_GLN)