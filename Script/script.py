#!/usr/bin/env python

import sys
import commandline as cl
import atom
import intcalc

if __name__ == "__main__":
    # check if the help is called
    if cl.check_help(sys.argv) == True: sys.exit()
    # check for potential errors in the command line
    if cl.check_command(sys.argv) == 1: sys.exit("A pdb file or arguments are missing")
    elif cl.check_command(sys.argv) == 2: sys.exit("A pdb file is needed")
    elif cl.check_command(sys.argv) == 3: sys.exit("No interactions calculation are given")

    # interaction using distance threshold and the default values
    interac_with_range = ["--hphb", "--inic", "--arar", "--arsu", "--capi"]
    interac_range = [5, 6, 4.5, 5.3, 6]
    # run all the calculation arguments of the command line
    for arg in range(2,len(sys.argv)):
        # read the pdb file
        with open(sys.argv[1], "r") as finput:
            # check if the calculation of the interaction need a distance threshold
            if sys.argv[arg][0:6] in interac_with_range:
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
                # calculation of ionic interactions
                sulphur = []
                # parsing of the pdb file
                for ligne in finput:
                    if ligne[0:6].strip() == "ATOM" and ligne[17:20].strip() == "CYS" and ligne[76:78].strip() == "S":
                        sulphur.append(atom.Atom(ligne[76:78].strip(), ligne[22:26].strip(), ligne[17:20].strip(),
                                                 ligne[21:22].strip(), ligne[30:38].strip(), ligne[38:46].strip(),
                                                 ligne[46:54].strip()))
                print("{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}".format("Position", "Residue", "Chain", "Position",
                                              "Residue", "Chain", "Distance"))
                # check the pair of sulphur
                for i, elem1 in enumerate(sulphur):
                    for elem2 in sulphur[i+1:]:
                        dist = intcalc.calc_range3D(elem1, elem2)
                        if dist <= 2.2:
                            print("{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10.2f}".format(elem1.position, elem1.residue,
                                                                                      elem1.chain, elem2.position,
                                                                                      elem2.residue, elem2.chain, dist))
                print("\n")

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