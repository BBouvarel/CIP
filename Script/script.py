#!/usr/bin/env python

import sys
import commandline as cl
import intcalc as ic

if __name__ == "__main__":
    # check if the help is called
    if cl.check_help(sys.argv) == True: sys.exit()
    # check for potential errors in the command line
    if cl.check_command(sys.argv) == 1: sys.exit("A pdb file or arguments are missing")
    elif cl.check_command(sys.argv) == 2: sys.exit("A pdb file is needed")
    elif cl.check_command(sys.argv) == 3: sys.exit("No interactions calculation are given")

    # interaction using distance threshold and the default values
    INTERAC_WITH_RANGE = ["--hphb", "--inic", "--arar", "--arsu", "--capi"]
    INTERAC_RANGE = [5, 6, 4.5, 5.3, 6]

    # run all the calculation arguments of the command line
    for arg in range(2, len(sys.argv)):
        # read the pdb file
        #with open(sys.argv[1], "r") as finput:
        # check if the calculation of the interaction need a distance threshold
        if sys.argv[arg][0:6] in INTERAC_WITH_RANGE:
            def_range = cl.set_val_default(sys.argv[arg], sys.argv[arg][0:6],
                                           INTERAC_RANGE[INTERAC_WITH_RANGE.index(sys.argv[arg][0:6])])
            if def_range == -1:
                sys.exit()

        if sys.argv[arg][0:6] == "--hphb":
            # calculation of hydrophobic interaction
            hphb = ic.parsing(sys.argv[1], ["CB", "CD", "CE", "CG", "CH", "CZ", "NE", "OH", "SD"],
                              ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"], [])
            #for i in range(0,len(hphb)): print(hphb[i].position)
            print("{:>45} {}A ".format("Hydrophobic interactions", def_range))
            print("{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}".format("Position",
                                                                      "Residue", "Chain", "Position",
                                                                      "Residue", "Chain", "Distance"))
            pos_prec = []
            # check the pair of sulphur
            for i, elem1 in enumerate(hphb):
                for elem2 in hphb[i + 1:]:
                    dist = ic.calc_range3D(elem1, elem2)
                    if dist <= def_range and ([elem1.position, elem2.position] not in pos_prec)\
                            and (elem1.position != elem2.position or elem1.chain != elem2.chain):
                        print("{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10.2f}".format(elem1.position, elem1.residue,
                                                                                     elem1.chain, elem2.position,
                                                                                     elem2.residue, elem2.chain, dist))
                        pos_prec.append([elem1.position, elem2.position])
            print("\n")

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
            # calculation of disulphide bridges
            sulphur = ic.parsing(sys.argv[1], [], ["CYS"], ["S"])
            print("{:^70}".format("Disulphide bridges"))
            print("{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}".format("Position", "Residue", "Chain", "Position",
                                                                      "Residue", "Chain", "Distance"))
            # check the pair of sulphur
            for i, elem1 in enumerate(sulphur):
                for elem2 in sulphur[i+1:]:
                    dist = ic.calc_range3D(elem1, elem2)
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

        else:
            sys.exit("L'un des arguments n'est pas reconnu par le programme")