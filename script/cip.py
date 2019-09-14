#!/usr/bin/env python

__author__ = "Bertrand Bouvarel"
__date__ = "2019/09"

import sys
import packages.commandline as cl
import packages.intcalc as ic

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
    INTERAC_RANGE = [5, 6, [4.5, 7], 5.3, 6]
    # run all the calculation arguments of the command line

    with open("../results/"+sys.argv[1][-8:-4]+"_res.txt", "w") as fout:
        fout.write("Results of the intra-protein interaction calculation:")
        # Create an empty result file

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
            #print("{:>45} {}A ".format("Hydrophobic interactions", def_range))
            ic.print_header(sys.argv[1][-8:-4], "Hydrophobic interactions", def_range)

            hphb = ic.parsing(sys.argv[1], ["CB", "CD", "CD1", "CD2", "CE", "CE1", "CE2", "CE3",
                                            "CG", "CG1", "CG2", "CH2", "CZ", "CZ2", "CZ3", "NE1", "SD"],
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
                        ic.print_pos_res_ch_dis(sys.argv[1][-8:-4], elem1.position, elem1.residue,
                                                elem1.chain, elem2.position,
                                                elem2.residue, elem2.chain, dist)
                        # print the results
                        pos_prev.append([elem1.position, elem2.position])
            print("\n")

        elif sys.argv[arg][0:6] == "--inic":
            # calculation of ionic interaction
            #print("{:>45} {}A ".format("Ionic interactions", def_range))
            ic.print_header(sys.argv[1][-8:-4], "Ionic interactions", def_range)

            inic = ic.parsing(sys.argv[1], ["ND1", "NH2", "NZ", "OD2", "OE2"],
                              ["ARG", "LYS", "HIS", "ASP", "GLU"])


            # ND1, NH2, NZ, OE2, OD2
            # "ND1", "ND2", "NH1", "NH2", "NZ", "OD1", "OD2", "OE1", "OE2"
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
                            ic.print_pos_res_ch_dis(sys.argv[1][-8:-4], elem1.position, elem1.residue,
                                                    elem1.chain, elem2.position,
                                                    elem2.residue, elem2.chain, dist)
                            pos_prev.append([elem1.position, elem2.position])
            print("\n")

        elif sys.argv[arg][0:6] == "--arar":
            # calculation of aromatic-aromatic interaction
            #print("{:>45} {}-{}A ".format("Aromatic-Aromatic interactions", def_range[0], def_range[1]))
            ic.print_header(sys.argv[1][-8:-4], "Aromatic-aromatic interactions", def_range)
            arar = ic.parsing(sys.argv[1], ["CD1", "CD2", "CE1", "CE2", "CG", "CZ"], ["PHE", "TYR"])
            arar = arar + ic.parsing(sys.argv[1], ["CD2", "CE2", "CE3", "CH2", "CZ2", "CZ3"], ["TRP"])

            pos_prev = []
            i = 0
            while i < len(arar):
                aro1 = ic.calc_centroid(arar[i:(i + 6)])
                j = i + 6
                while j < len(arar):
                    aro2 = ic.calc_centroid(arar[j:(j+6)])
                    dist = ic.calc_range(aro1, aro2)
                    #if arar[i].position == 29 and arar[j].position == 44: print(dist, aro1.residue, aro2.residue)
                    if def_range[0] <= dist <= def_range[1] and \
                            ([aro1.position, aro2.position] not in pos_prev) and \
                            (aro1.position != aro2.position or
                             aro1.chain != aro2.chain):
                        ic.print_pos_res_ch_dis(sys.argv[1][-8:-4], aro1.position, aro1.residue,
                                                aro1.chain, aro2.position,
                                                aro2.residue, aro2.chain, dist)
                        pos_prev.append([aro1.position, aro2.position])
                    j += 6
                i += 6
            print("\n")

        elif sys.argv[arg][0:6] == "--arsu":
            # calculation of aromatic-sulphur interaction
            #print("{:>45} {}A ".format("Aromatic-Sulphure interactions", def_range))
            ic.print_header(sys.argv[1][-8:-4], "Aromatic-sulphure interactions", def_range)
            aro_all = ic.parsing(sys.argv[1], ["CD1", "CD2", "CE1", "CE2", "CG", "CZ"], ["PHE", "TYR"])
            aro_all = aro_all + ic.parsing(sys.argv[1], ["CD2", "CE2", "CE3", "CH2", "CZ2", "CZ3"], ["TRP"])
            sul_all = ic.parsing(sys.argv[1], ["SD", "SG"], ["CYS", "MET"])
            pos_prev = []
            i = 0
            while i < len(aro_all):
                aro = ic.calc_centroid(aro_all[i:(i + 6)])
                for sul in sul_all:
                    dist = ic.calc_range(aro, sul)
                    if dist <= def_range and \
                            ([aro.position, sul.position] not in pos_prev) and \
                            (aro.position != sul.position or
                             aro.chain != sul.chain):
                        ic.print_pos_res_ch_dis(sys.argv[1][-8:-4], aro.position, aro.residue,
                                                aro.chain, sul.position,
                                                sul.residue, sul.chain, dist)
                        pos_prev.append([aro.position, sul.position])
                i += 6
            print("\n")


        elif sys.argv[arg][0:6] == "--capi":
            # calculation of Cation-pi interaction
            #print("{:>45} {}A ".format("Cation-pi interactions", def_range))
            ic.print_header(sys.argv[1][-8:-4], "Cation-pi interactions", def_range)
            #aro_all = ic.parsing(sys.argv[1], ["CB", "CD1", "CD2", "CE1", "CE2", "CE3", "CG",
             #                                  "CH2", "CZ", "CZ2", "CZ3", "NE1", "OH"],
             #                    ["PHE", "TRP", "TYR"])
            aro_all = ic.parsing(sys.argv[1], ["CD1", "CD2", "CE1", "CE2", "CG", "CZ"], ["PHE", "TYR"])
            aro_all = aro_all + ic.parsing(sys.argv[1], ["CD2", "CE2", "CE3", "CH2", "CZ2", "CZ3"], ["TRP"])
            cation_all = ic.parsing(sys.argv[1], ["NH2", "NZ"], ["ARG", "LYS"])
            # cation_all = ic.parsing(sys.argv[1], ["CB", "CD", "CE", "CG", "CZ", "NE", "NH1",
            #                                                   "NH2", "NZ"], ["ARG", "LYS"])

            pos_prev = []
            i = 0
            while i < len(aro_all):
                aro = ic.calc_centroid(aro_all[i:(i + 6)])
                for ato_cat in cation_all:
                    dist = ic.calc_range(aro, ato_cat)
                    if dist <= def_range and \
                            ([aro.position, ato_cat.position] not in pos_prev) and \
                            (aro.position != ato_cat.position or
                             aro.chain != ato_cat.chain):
                        ic.print_pos_res_ch_dis(sys.argv[1][-8:-4], aro.position, aro.residue,
                                                aro.chain, ato_cat.position,
                                                ato_cat.residue, ato_cat.chain, dist)
                        pos_prev.append([aro.position, ato_cat.position])
                i += 6
            print("\n")

        elif sys.argv[arg][0:6] == "--disu":
            # calculation of disulphide bridges
            #print("{:^70}".format("Disulphide bridges 2.2A"))
            ic.print_header(sys.argv[1][-8:-4], "Disulphide bridges", 2.2)

            sulphur = ic.parsing(sys.argv[1], ["SG"], ["CYS"])
            for i, elem1 in enumerate(sulphur):
                # check the pair of sulphur
                for elem2 in sulphur[i + 1:]:
                    dist = ic.calc_range(elem1, elem2)
                    if dist <= 2.2:
                        ic.print_pos_res_ch_dis(sys.argv[1][-8:-4], elem1.position, elem1.residue,
                                                elem1.chain, elem2.position,
                                                elem2.residue, elem2.chain, dist)
            print("\n")

        elif sys.argv[arg][0:6] == "--mmhb":
            # calculation of hydrogen bonds main-main
            # print("{:^90}\n{:^40}{:^40}".format("Hydrogen bonds main-main", "Donors", "Acceptors"))
            ic.print_hydrogen_header(sys.argv[1][-8:-4], "Hydrogen bonds main-main", [3.5, 4])

            pos_prev = []
            mmhb = ic.parsing(sys.argv[1], ["N", "O", "OXT"], [])
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
                        ic.print_hydrogen_res(sys.argv[1][-8:-4], donor.position, donor.residue,
                                              donor.chain, donor.name,
                                              acceptor.position, acceptor.residue,
                                              acceptor.chain, acceptor.name, dist)
                        pos_prev.append([donor.position, acceptor.position])
            print("\n")

        elif sys.argv[arg][0:6] == "--mshb":
            # calculation of hydrogen bonds main-side
            #print("{:^90}\n{:^40}{:^40}".format("Hydrogen bonds main-side", "Donors", "Acceptors"))
            ic.print_hydrogen_header(sys.argv[1][-8:-4], "Hydrogen bonds main-side", [3.5, 4])

            main = ["N", "O", "OXT"]
            #side = ["OD", "OE", "OG", "OH", "ND", "NE", "NH", "NZ", "SD", "SG"]
            side_all = ["ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ", "OD1", "OD2", "OE1",
                        "OE2", "OG", "OG1", "OH", "SD", "SG"]

            side_don = {"ARG": ["NE", "NH1", "NH2"], "ASN": ["ND2"], "CYS": ["SG"],
                        "GLN": ["NE2"], "HIS": ["ND1", "NE2"], "LYS": ["NZ"],
                        "SER": ["OG"], "THR": ["OG1"], "TRP": ["NE1"], "TYR": ["OH"]}
            side_acc = {"ASN": ["OD1"], "ASP": ["OD1", "OD2"], "GLN": ["OE1"],
                        "GLU": ["OE1", "OE2"], "HIS": ["ND1", "NE2"], "MET": ["SD"], "SER": ["OG"],
                        "THR": ["OG1"], "TYR": ["OH"]}
            pos_prev = []
            donor = None
            acceptor = None
            prev_do_ac = (None, None)
            mshb = ic.parsing(sys.argv[1], main + side_all, [])

            for i, elem1 in enumerate(mshb):
                for elem2 in mshb[i + 1:]:
                    def_range = 3.5
                    #print(elem1.name, elem2.name)
                    if elem1.name in main and elem2.name in side_all:
                        if elem1.name == "N" and\
                                (elem2.residue in side_acc.keys() and elem2.name in side_acc[elem2.residue]):
                            donor = elem1
                            acceptor = elem2
                        elif (elem1.name == "O" or elem1.name == "OXT") and\
                                (elem2.residue in side_don.keys() and elem2.name in side_don[elem2.residue]):
                            donor = elem2
                            acceptor = elem1
                    elif elem1.name in side_all and elem2.name in main:
                        if elem2.name == "N" and \
                                (elem1.residue in side_acc.keys() and elem1.name in side_acc[elem1.residue]):
                            donor = elem2
                            acceptor = elem1
                        elif (elem2.name == "O" or elem2.name == "OXT") and \
                                (elem1.residue in side_don.keys() and elem1.name in side_don[elem1.residue]):
                            donor = elem1
                            acceptor = elem2
                    else:
                        continue

                    if donor is not None and acceptor is not None and (donor, acceptor) != prev_do_ac:
                        if donor.residue == "PRO" and donor.name == "N":
                            continue
                        if donor.name == "SG" or acceptor.name == "SD":
                            def_range = 4
                        dist = ic.calc_range(donor, acceptor)
                        if dist <= def_range and \
                                ([donor.position, acceptor.position] not in pos_prev or donor.name != acceptor.name) and \
                                (donor.position != acceptor.position or
                                 donor.chain != acceptor.chain or
                                 donor.name != acceptor.name) and \
                                abs(donor.position - acceptor.position) >= 2:

                            ic.print_hydrogen_res(sys.argv[1][-8:-4], donor.position, donor.residue,
                                                  donor.chain, donor.name,
                                                  acceptor.position, acceptor.residue,
                                                  acceptor.chain, acceptor.name, dist)
                            pos_prev.append([donor.position, acceptor.position])
                        prev_do_ac = (donor, acceptor)
            print("\n")

        elif sys.argv[arg][0:6] == "--sshb":
            # calculation of hydrogen bonds side-side
            #print("{:^90}\n{:^40}{:^40}".format("Hydrogen bonds side-side", "Donors", "Acceptors"))
            ic.print_hydrogen_header(sys.argv[1][-8:-4], "Hydrogen bonds side-side", [3.5, 4])

            side_all = ["ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ", "OD1", "OD2", "OE1",
                        "OE2", "OG", "OG1", "OH", "SD", "SG"]
            poss_don = {"ARG": ["NE", "NH1", "NH2"], "ASN": ["ND2"], "CYS": ["SG"],
                        "GLN": ["NE2"], "HIS": ["ND1", "NE2"], "LYS": ["NZ"],
                        "SER": ["OG"], "THR": ["OG1"], "TRP": ["NE1"], "TYR": ["OH"]}
            poss_acc = {"ASN": ["OD1"], "ASP": ["OD1", "OD2"], "GLN": ["OE1"],
                        "GLU": ["OE1", "OE2"], "HIS": ["ND1", "NE2"], "MET": ["SD"], "SER": ["OG"],
                        "THR": ["OG1"], "TYR": ["OH"]}
            donor = None
            acceptor = None
            pos_prev = []
            prev_do_ac = (None, None)
            sshb = ic.parsing(sys.argv[1], side_all, [])
            for i, elem1 in enumerate(sshb):
                for elem2 in sshb:
                    def_range = 3.5
                    if (elem1.residue in poss_don.keys() and elem1.name in poss_don[elem1.residue]) and\
                            (elem2.residue in poss_acc.keys() and elem2.name in poss_acc[elem2.residue]):
                        donor = elem1
                        acceptor = elem2
                    elif (elem1.residue in poss_acc.keys() and elem1.name in poss_acc[elem1.residue]) and\
                            (elem2.residue in poss_don.keys() and elem2.name in poss_don[elem2.residue]):
                        donor = elem2
                        acceptor = elem1
                    else:
                        continue
                    if donor.name == "SG" or acceptor.name == "SD":
                        def_range = 4
                    dist = ic.calc_range(donor, acceptor)
                    if dist <= def_range and \
                            ([donor.position, acceptor.position] not in pos_prev) and \
                            (donor.position != acceptor.position or
                             donor.chain != acceptor.chain or
                             donor.name != acceptor.name) and \
                            abs(donor.position - acceptor.position) >= 2:
                        ic.print_hydrogen_res(sys.argv[1][-8:-4], donor.position, donor.residue,
                                              donor.chain, donor.name,
                                              acceptor.position, acceptor.residue,
                                              acceptor.chain, acceptor.name, dist)
                        pos_prev.append([donor.position, acceptor.position])
            print("\n")

        else:
            sys.exit("One of the arguments is not recognized by the program")
