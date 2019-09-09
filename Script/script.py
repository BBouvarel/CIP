#!/usr/bin/env python

import sys
import command_line as cl
import math
import Atome


# faire un if argv = 2 (script.py --help),alors afficher help des commande
# else le reste

if __name__ == "__main__":
    if cl.check_help(sys.argv) == True: sys.exit()
    if cl.check_command(sys.argv) == 1: sys.exit("Il vous manque un fichier pdb et un ou plusieurs arguments")
    elif cl.check_command(sys.argv) == 2: sys.exit("Il est necessaire de transmettre un fichier pdb")
    elif cl.check_command(sys.argv) == 3: sys.exit("Aucun calcul d'interaction specifie")

    interac_with_range = ["--hphb", "--inic", "--arar", "--arsu", "--capi"]
    interac_range = [5, 6, 4.5, 5.3, 6]
    for arg in range(2,len(sys.argv)):
        # Parcour de l'ensemble des arguments transmis dans la commande

        if sys.argv[arg][0:6] in interac_with_range:
            # print(interac_range.index(sys.argv[arg][0:6]), "aaa")
            def_range = cl.set_val_default(sys.argv[arg], sys.argv[arg][0:6],
                                           interac_range[interac_with_range.index(sys.argv[arg][0:6])])
            if def_range == -1:
                sys.exit()

        if sys.argv[arg][0:6] == "--hphb":
            print(def_range)
        # Lancer calcul interac hydrophobicite

        elif sys.argv[arg][0:6] == "--inic":
            print(def_range)
            # Lancer calcul interac ionic

        elif sys.argv[arg][0:6] == "--arar":
            print(def_range)
            # Lancer calcul interac aromatic-aromatic

        elif sys.argv[arg][0:6] == "--arsu":
            print(def_range)
            # Lancer calcul interac aromatic-sulphur

        elif sys.argv[arg][0:6] == "--capi":
            print(def_range)
            # Lancer calcul interac cation-pi

        elif sys.argv[arg][0:6] == "--disu":
            # Lancer calcul disulphide bridges
            print("disul")

        elif sys.argv[arg][0:6] == "--mmhb":
            # Lancer calcul interac hydrogen main chain - main chain
            print("mmhb")

        elif sys.argv[arg][0:6] == "--mshb":
            # Lancer calcul interac hydrogen main chain - side chain
            print("mshb")

        elif sys.argv[arg][0:6] == "--sshb":
            # Lancer calcul interac hydrogen side chain - side chain
            print("sshb")

        else :
            sys.exit("L'un des arguments n'est pas reconnu par le programme")