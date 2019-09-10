#!/usr/bin/env python

import sys
import commandline as cl
import intcalc as ic

if __name__ == "__main__":
    # check if the help is called
    if cl.check_help(sys.argv):
        sys.exit()
    # check for potential errors in the command line
    if cl.check_command(sys.argv) == 1:
        sys.exit("A pdb file or arguments are missing")
    elif cl.check_command(sys.argv) == 2:
        sys.exit("A pdb file is needed")
    elif cl.check_command(sys.argv) == 3:
        sys.exit("No interactions calculation are given")

    # interaction using distance threshold and the default values
    INTERAC_WITH_RANGE = ["--hphb", "--inic", "--arar", "--arsu", "--capi"]
    INTERAC_RANGE = [5, 6, 4.5, 5.3, 6]
    # run all the calculation arguments of the command line
    for arg in range(2, len(sys.argv)):
        # check if the calculation of the interaction need a distance threshold
        if sys.argv[arg][0:6] in INTERAC_WITH_RANGE:
            def_range = cl.set_val_default(sys.argv[arg], sys.argv[arg][0:6],
                                           INTERAC_RANGE[INTERAC_WITH_RANGE.index(sys.argv[arg][0:6])])
            if def_range == -1:
                sys.exit()

        if sys.argv[arg][0:6] == "--hphb":
            # calculation of hydrophobic interaction
            print("{:>45} {}A ".format("Hydrophobic interactions", def_range))
            ic.print_results()

            hphb = ic.parsing(sys.argv[1], ["CB", "CD", "CE", "CG", "CH", "CZ", "NE", "OH", "SD"],
                              ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"])
            pos_prec = []
            for i, elem1 in enumerate(hphb):
                # check the pair of elements
                for elem2 in hphb[i + 1:]:
                    dist = ic.calc_range(elem1, elem2)
                    if dist <= def_range and\
                            ([elem1.get_position(), elem2.get_position()] not in pos_prec) and\
                            (elem1.get_position() != elem2.get_position() or
                             elem1.get_chain() != elem2.get_chain()):
                        #ic.print_pos_res_ch_dis(elem1.get_position(), elem1.get_residue(),
                         #                        elem1.get_chain(), elem2.get_position(),
                          #                       elem2.get_residue(), elem2.get_chain(), dist)
                        print("{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}"
                              "{:^10.2f}".format(elem1.get_position(),
                                                 elem1.get_residue(),
                                                 elem1.get_chain(),
                                                 elem2.get_position(),
                                                 elem2.get_residue(),
                                                 elem2.get_chain(),
                                                 dist))
                        pos_prec.append([elem1.get_position(), elem2.get_position()])
            print("\n")

        elif sys.argv[arg][0:6] == "--inic":
            # calculation of ionic interaction
            print("{:>45} {}A ".format("Ionic interactions", def_range))
            ic.print_results()

            inic = ic.parsing(sys.argv[1], ["ND", "NE", "NH", "NZ", "OD", "OE"],
                              ["ARG", "LYS", "HIS", "ASP", "GLU"])
            pos_prec = []
            pos_res = ["ARG", "LYS", "HIS"]
            neg_res = ["ASP", "GLU"]
            for i, elem1 in enumerate(inic):
                # check the pair of elements
                for elem2 in inic[i + 1:]:
                    dist = ic.calc_range(elem1, elem2)
                    if dist <= def_range and\
                            ([elem1.get_position(), elem2.get_position()] not in pos_prec) and\
                            (elem1.get_position() != elem2.get_position() or
                             elem1.get_chain() != elem2.get_chain()):
                        if (elem1.get_residue() in pos_res and elem2.get_residue() in neg_res) or\
                                (elem1.get_residue() in neg_res and elem2.get_residue() in pos_res):
                            print("{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}"
                                  "{:^10.2f}".format(elem1.get_position(),
                                                     elem1.get_residue(),
                                                     elem1.get_chain(),
                                                     elem2.get_position(),
                                                     elem2.get_residue(),
                                                     elem2.get_chain(),
                                                     dist))
                            pos_prec.append([elem1.get_position(), elem2.get_position()])
            print("\n")

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
            print("{:^70}".format("Disulphide bridges 2.2A"))
            ic.print_results()

            sulphur = ic.parsing(sys.argv[1], ["SG"], ["CYS"])
            for i, elem1 in enumerate(sulphur):
                # check the pair of sulphur
                for elem2 in sulphur[i+1:]:
                    dist = ic.calc_range(elem1, elem2)
                    if dist <= 2.2:
                        print("{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}"
                              "{:^10.2f}".format(elem1.get_position(),
                                                 elem1.get_residue(),
                                                 elem1.get_chain(),
                                                 elem2.get_position(),
                                                 elem2.get_residue(),
                                                 elem2.get_chain(),
                                                 dist))
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
