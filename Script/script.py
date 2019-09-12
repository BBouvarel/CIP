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

        # LANE DOUBLE BOUCLE PUIS FI CHECK LA COMMANDE QUI PERMET DE LANCER LA FONC DE L'INTERAC QUI CORRESPOND AVEC LES PRINT DANS LA FONCTION
        if sys.argv[arg][0:6] == "--hphb":
            # calculation of hydrophobic interaction
            print("{:>45} {}A ".format("Hydrophobic interactions", def_range))
            ic.print_header()

            hphb = ic.parsing(sys.argv[1], ["CB", "CD", "CE", "CG", "CH", "CZ", "NE", "OH", "SD"],
                              ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"])
            pos_prev = []
            for i, elem1 in enumerate(hphb):
                for elem2 in hphb[i + 1:]:
                    # check the pairs of elements
                    dist = ic.calc_range(elem1, elem2)
                    if dist <= def_range and \
                            ([elem1.position, elem2.position] not in pos_prev) and \
                            (elem1.position != elem2.position or
                             elem1.chain != elem2.chain):
                        ic.print_pos_res_ch_dis(elem1.position, elem1.residue,
                                                elem1.chain, elem2.position,
                                                elem2.residue, elem2.chain, dist)
                        # print the results
                        pos_prev.append([elem1.position, elem2.position])
            print("\n")

        elif sys.argv[arg][0:6] == "--inic":
            # calculation of ionic interaction
            print("{:>45} {}A ".format("Ionic interactions", def_range))
            ic.print_header()

            inic = ic.parsing(sys.argv[1], ["ND", "NE", "NH", "NZ", "OD", "OE"],
                              ["ARG", "LYS", "HIS", "ASP", "GLU"])
            pos_prev = []
            pos_res = ["ARG", "LYS", "HIS"]
            neg_res = ["ASP", "GLU"]
            for i, elem1 in enumerate(inic):
                # check the pair of elements
                for elem2 in inic[i + 1:]:
                    dist = ic.calc_range(elem1, elem2)
                    if dist <= def_range and \
                            ([elem1.position, elem2.position] not in pos_prev) and \
                            (elem1.position != elem2.position or
                             elem1.chain != elem2.chain):
                        if (elem1.residue in pos_res and elem2.residue in neg_res) or \
                                (elem1.residue in neg_res and elem2.residue in pos_res):
                            # binding of a positive res with a negative res only
                            ic.print_pos_res_ch_dis(elem1.position, elem1.residue,
                                                    elem1.chain, elem2.position,
                                                    elem2.residue, elem2.chain, dist)
                            pos_prev.append([elem1.position, elem2.position])
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
            ic.print_header()

            sulphur = ic.parsing(sys.argv[1], ["SG"], ["CYS"])
            for i, elem1 in enumerate(sulphur):
                # check the pair of sulphur
                for elem2 in sulphur[i + 1:]:
                    dist = ic.calc_range(elem1, elem2)
                    if dist <= 2.2:
                        ic.print_pos_res_ch_dis(elem1.position, elem1.residue,
                                                elem1.chain, elem2.position,
                                                elem2.residue, elem2.chain, dist)
            print("\n")

        elif sys.argv[arg][0:6] == "--mmhb":
            # calculation of hydrogen bonds main-main
            print("{:^90}\n{:^40}{:^40}".format("Hydrogen bonds main-main", "Donnor", "Acceptor"))
            ic.print_hydrogen_header()

            pos_prev = []
            mmhb = ic.parsing(sys.argv[1], ["N", "O"], [])
            for i, elem1 in enumerate(mmhb):
                for elem2 in mmhb[i + 1:]:
                    if elem1.name == "N" and elem2.name == "O":
                        donor = elem1
                        acceptor = elem2
                    else:
                        donor = elem2
                        acceptor = elem1
                    dist = ic.calc_range(donor, acceptor)
                    if dist <= 3.5 and \
                            ([donor.position, acceptor.position] not in pos_prev) and \
                            (donor.position != acceptor.position or
                             donor.chain != acceptor.chain) and \
                            abs(donor.position - acceptor.position) >= 2 and \
                            donor.residue != "PRO":
                        ic.print_hydrogen_res(donor.position, donor.residue,
                                              donor.chain, donor.name,
                                              acceptor.position, acceptor.residue,
                                              acceptor.chain, acceptor.name, dist)
                        pos_prev.append([donor.position, acceptor.position])
            print("\n")

        elif sys.argv[arg][0:6] == "--mshb":
            # calculation of hydrogen bonds main-side
            print("{:^90}\n{:^40}{:^40}".format("Hydrogen bonds side-main", "Donor", "Acceptor"))
            ic.print_hydrogen_header()

            main = ["N", "O", "OXT"]
            side = ["OD", "OE", "OG", "OH", "ND", "NE", "NH", "NZ", "SG"]
            pos_prev = []
            donor = None
            acceptor = None
            do_ac = [(None, None)]
            mshb = ic.parsing(sys.argv[1], main + side, [])
            print(len(mshb))

            for i, elem1 in enumerate(mshb):
                for elem2 in mshb[i + 1:]:
                    def_range = 3.5

                    if (elem1.name in main and elem2.name in side):
                        if elem1.name == "N":
                            donor = elem1
                            acceptor = elem2
                        else:
                            donor = elem2
                            acceptor = elem1
                        if elem2.name == "SG":
                            def_range = 4
                    elif (elem1.name in side and elem2.name in main):
                        if elem2.name == "N":
                            donor = elem2
                            acceptor = elem1
                        else:
                            donor = elem1
                            acceptor = elem2
                        if elem1.name == "SG":
                            def_range = 4
                    if donor is not None and acceptor is not None and (donor, acceptor) not in do_ac:
                        dist = ic.calc_range(donor, acceptor)
                        if dist <= def_range and \
                                ([donor.position, acceptor.position] not in pos_prev or donor.name != acceptor.name) and \
                                (donor.position != acceptor.position or
                                 donor.chain != acceptor.chain or
                                 donor.name != acceptor.name) and \
                                abs(donor.position - acceptor.position) >= 2:
                            ic.print_hydrogen_res(donor.position, donor.residue,
                                                  donor.chain, donor.name,
                                                  acceptor.position, acceptor.residue,
                                                  acceptor.chain, acceptor.name, dist)
                            pos_prev.append([donor.position, acceptor.position])
                        do_ac.append((donor, acceptor))
            print("\n")

        elif sys.argv[arg][0:6] == "--sshb":
            # calculation of hydrogen bonds side-side
            print("{:^90}\n{:^40}{:^40}".format("Hydrogen bonds side-side", "Donnor", "Acceptor"))
            ic.print_hydrogen_header()

            pos_prev = []
            mmhb = ic.parsing(sys.argv[1], ["N", "O"], [])
            for i, elem1 in enumerate(mmhb):
                for elem2 in mmhb[i + 1:]:
                    if elem1.name == "N" and elem2.name == "O":
                        donor = elem1
                        acceptor = elem2
                    elif elem1.name == "O" and elem2.name == "N":
                        donor = elem2
                        acceptor = elem1
                    dist = ic.calc_range(donor, acceptor)
                    if dist <= 3.5 and \
                            ([donor.position, acceptor.position] not in pos_prev) and \
                            (donor.position != acceptor.position or
                             donor.chain != acceptor.chain) and \
                            abs(donor.position - acceptor.position) >= 2 and \
                            donor.residue != "PRO":
                        ic.print_hydrogen_res(donor.position, donor.residue,
                                              donor.chain, donor.name,
                                              acceptor.position, acceptor.residue,
                                              acceptor.chain, acceptor.name, dist)
                        pos_prev.append([donor.position, acceptor.position])
            print("\n")

        else:
            sys.exit("L'un des arguments n'est pas reconnu par le programme")
