#!/usr/bin/env python

"""
package containing all functions, calculations, display
and writing of protein interactions of a PDB file
"""

__author__ = "Bertrand BOUVAREL"
__date__ = "2019/09"

import math
import statistics as stat
import packages.atom as atom


def add_to_list(line, list):
    """
    function used to add an atom to the list of those meeting specific criteria,
    used to clarify the code

    :param line: line of the pdb file containing the element of interest
    :param list: list containing the objects of class atom
    :return: list of objects of class atom with the new element added
    """
    list.append(atom.Atom(line[12:16].strip(), line[22:26].strip(),
                          line[17:20].strip(), line[21:22].strip(),
                          line[30:38].strip(), line[38:46].strip(),
                          line[46:54].strip()))
    return list


def parsing(pdb, name, res):
    """
    function allowing the parsing of the pdb file following certain conditions

    :param pdb: the pdb file given by the user
    :param name: atom name
    :param res: list of residues containing the element
    :return: list of objects of class atom that respect the conditions
    """
    elements = []
    with open(pdb, "r") as finput:
        for line in finput:
            if line[0:6].strip() == "ATOM":
                # if the line describe an atom
                if line[17:20].strip() in res and res:
                    # if the residue of the line contain this specific atom
                    if line[12:16].strip() in name and name:
                        # if the atom of the line is the desired atom
                        elements = add_to_list(line, elements)
                elif not res:
                    # if there is no residue criterion
                    if line[12:16].strip() in name and name:
                        elements = add_to_list(line, elements)
    return elements


def calc_range(elem1, elem2):
    """
    function allowing the calculation of the distance of two atoms

    :param elem1: object of the class atom
    :param elem2: object of the class atom
    :return: distance between two atoms
    """
    dist = math.sqrt((elem2.coor_x-elem1.coor_x)**2 +
                     (elem2.coor_y-elem1.coor_y)**2 +
                     (elem2.coor_z-elem1.coor_z)**2)
    return dist


def calc_centroid(atoms):
    """
    function calculating the centroid of a phenyl group of an aromatic residue

    :param atoms: a list of objects of the atom class corresponding to the six atoms of a phenyl
    :return: an object containing the information about the phenyl and his centroid
    """
    return atom.Atom("Phenyl", atoms[0].position, atoms[0].residue, atoms[0].chain,
                     stat.mean([atoms[0].coor_x, atoms[1].coor_x, atoms[2].coor_x,
                                atoms[3].coor_x, atoms[4].coor_x, atoms[5].coor_x]),
                     stat.mean([atoms[0].coor_y, atoms[1].coor_y, atoms[2].coor_y,
                                atoms[3].coor_y, atoms[4].coor_y, atoms[5].coor_y]),
                     stat.mean([atoms[0].coor_z, atoms[1].coor_z, atoms[2].coor_z,
                                atoms[3].coor_z, atoms[4].coor_z, atoms[5].coor_z]))


def check_criteria(intra_inter, elem1, elem2, pos_prev, dist, def_range):
    """
    function checking the criteria for hydrophobic, ionic, aromatic-aromatic,
    aromatic-sulphur and cation-pi interactions, used to clarify the code

    :param intra_inter: choice to calculate intra or inter protein interactions
    :param elem1: object of class atom
    :param elem2: object of class atom
    :param pos_prev: list of pairs of element already printed
    :param dist: distance between the two elements
    :param def_range: threshold distance
    :return: a boolean to know if the criteria are respected
    """
    if dist <= def_range and \
            ([elem1.position, elem2.position, elem1.chain, elem2.chain] not in pos_prev) and \
            elem1.check_non_id(elem2):
        if intra_inter == "--intra" and elem1.chain == elem2.chain:
            return True
        elif intra_inter == "--inter" and elem1.chain != elem2.chain:
            return True
        else:
            return False
    else:
        return False


def hydrophobic(arg, def_range):
    """
    function calculating the hydrophobic interactions

    :param arg: the list of arguments of the command line
    :param def_range: the threshold distance, criterion of the interaction
    """
    print_header(arg[1][-8:-4], "Hydrophobic interactions", def_range)
    hphb = parsing(arg[1],
                   ["CB", "CD", "CD1", "CD2", "CE", "CE1", "CE2", "CE3", "CG", "CG1", "CG2", "CH2",
                    "CZ", "CZ2", "CZ3", "SD"],
                   ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"])

    pos_prev = []
    for i, elem1 in enumerate(hphb):
        for elem2 in hphb[i + 1:]:
            # check the pairs of elements
            dist = calc_range(elem1, elem2)
            if check_criteria(arg[2], elem1, elem2, pos_prev, dist, def_range):
                print_pos_res_ch_dis(arg[1][-8:-4], elem1.position, elem1.residue, elem1.chain,
                                     elem2.position, elem2.residue, elem2.chain, dist)
                pos_prev.append([elem1.position, elem2.position, elem1.chain, elem2.chain])
    print("\n")


def ionic(arg, def_range):
    """
    function calculating the ionic interactions

    :param arg: the list of arguments of the command line
    :param def_range: the threshold distance, criterion of the interaction
    """
    print_header(arg[1][-8:-4], "Ionic interactions", def_range)
    inic = parsing(arg[1], ["ND1", "NE", "NH1", "NH2", "NZ", "NE2", "OD2", "OE2"],
                   ["ARG", "LYS", "HIS", "ASP", "GLU"])
    # NE and NH1 because of the displacement of the electron in the guanidium group of the arginine
    pos_prev = []
    pos_res = ["ARG", "LYS", "HIS"]
    neg_res = ["ASP", "GLU"]
    for i, elem1 in enumerate(inic):
        for elem2 in inic[i + 1:]:
            # check the pairs of elements
            dist = calc_range(elem1, elem2)
            if check_criteria(arg[2], elem1, elem2, pos_prev, dist, def_range):
                if (elem1.residue in pos_res and elem2.residue in neg_res) or \
                        (elem1.residue in neg_res and elem2.residue in pos_res):
                    # binding of a positive res with a negative res only
                    print_pos_res_ch_dis(arg[1][-8:-4], elem1.position, elem1.residue, elem1.chain,
                                         elem2.position, elem2.residue, elem2.chain, dist)
                    pos_prev.append([elem1.position, elem2.position, elem1.chain, elem2.chain])
    print("\n")


def aro_aro(arg, def_range):
    """
    function calculating the aromatic-aromatic interactions

    :param arg: the list of arguments of the command line
    :param def_range: the threshold distance, criterion of the interaction
    """
    print_header(arg[1][-8:-4], "Aromatic-aromatic interactions", def_range)
    arar = parsing(arg[1], ["CD1", "CD2", "CE1", "CE2", "CG", "CZ"], ["PHE", "TYR"])
    arar = arar + parsing(arg[1], ["CD2", "CE2", "CE3", "CH2", "CZ2", "CZ3"], ["TRP"])
    pos_prev = []
    i = 0
    while i < len(arar):
        aro1 = calc_centroid(arar[i:(i + 6)])
        j = i + 6
        while j < len(arar):
            # check the pairs of aromatic rings
            aro2 = calc_centroid(arar[j:(j + 6)])
            dist = calc_range(aro1, aro2)
            if check_criteria(arg[2], aro1, aro2, pos_prev, dist, def_range[1]):
                if def_range[0] <= dist:
                    print_pos_res_ch_dis(arg[1][-8:-4], aro1.position, aro1.residue, aro1.chain,
                                         aro2.position, aro2.residue, aro2.chain, dist)
                    pos_prev.append([aro1.position, aro2.position, aro1.chain, aro2.chain])
            j += 6
        i += 6
    print("\n")


def aro_sul(arg, def_range):
    """
    function calculating the aromatic-sulphur interactions

    :param arg: the list of arguments of the command line
    :param def_range: the threshold distance, criterion of the interaction
    """
    print_header(arg[1][-8:-4], "Aromatic-sulphur interactions", def_range)
    aro_all = parsing(arg[1], ["CD1", "CD2", "CE1", "CE2", "CG", "CZ"], ["PHE", "TYR"])
    aro_all = aro_all + parsing(arg[1], ["CD2", "CE2", "CE3", "CH2", "CZ2", "CZ3"], ["TRP"])
    sul_all = parsing(arg[1], ["SD", "SG"], ["CYS", "MET"])
    pos_prev = []
    i = 0
    while i < len(aro_all):
        aro = calc_centroid(aro_all[i:(i + 6)])
        for sul in sul_all:
            # check the aromatic-sulphur pairs
            dist = calc_range(aro, sul)
            if check_criteria(arg[2], aro, sul, pos_prev, dist, def_range):
                print_pos_res_ch_dis(arg[1][-8:-4], aro.position, aro.residue, aro.chain,
                                     sul.position, sul.residue, sul.chain, dist)
                pos_prev.append([aro.position, sul.position, aro.chain, sul.chain])
        i += 6
    print("\n")


def cation_pi(arg, def_range):
    """
    function calculating the cation-pi interactions

    :param arg: the list of arguments of the command line
    :param def_range: the threshold distance, criterion of the interaction
    """
    print_header(arg[1][-8:-4], "Cation-pi interactions", def_range)
    aro_all = parsing(arg[1], ["CD1", "CD2", "CE1", "CE2", "CG", "CZ"], ["PHE", "TYR"])
    aro_all = aro_all + parsing(arg[1], ["CD2", "CE2", "CE3", "CH2", "CZ2", "CZ3"], ["TRP"])
    cation_all = parsing(arg[1], ["NH2", "NZ"], ["ARG", "LYS"])
    pos_prev = []
    i = 0
    while i < len(aro_all):
        aro = calc_centroid(aro_all[i:(i + 6)])
        for ato_cat in cation_all:
            # check the positive atom-aromatic pairs
            dist = calc_range(aro, ato_cat)
            if check_criteria(arg[2], aro, ato_cat, pos_prev, dist, def_range):
                print_pos_res_ch_dis(arg[1][-8:-4], aro.position, aro.residue, aro.chain,
                                     ato_cat.position, ato_cat.residue, ato_cat.chain, dist)
                pos_prev.append([aro.position, ato_cat.position, aro.chain, ato_cat.chain])
        i += 6
    print("\n")


def disulphide(arg):
    """
    function calculating the disulphide bridges

    :param arg: the list of arguments of the command line
    """
    print_header(arg[1][-8:-4], "Disulphide bridges", 2.2)
    sulphur = parsing(arg[1], ["SG"], ["CYS"])
    for i, elem1 in enumerate(sulphur):
        # check the pair of sulphur
        for elem2 in sulphur[i + 1:]:
            dist = calc_range(elem1, elem2)
            if dist <= 2.2:
                if arg[2] == "--intra" and elem1.chain == elem2.chain:
                    print_pos_res_ch_dis(arg[1][-8:-4], elem1.position, elem1.residue, elem1.chain,
                                         elem2.position, elem2.residue, elem2.chain, dist)
                elif arg[2] == "--inter" and elem1.chain != elem2.chain:
                    print_pos_res_ch_dis(arg[1][-8:-4], elem1.position, elem1.residue, elem1.chain,
                                         elem2.position, elem2.residue, elem2.chain, dist)
                else:
                    continue
    print("\n")


def mm_hbond(arg):
    """
    function calculating the main chain-main chain hydrogen bonds

    :param arg: the list of arguments of the command line
    """
    print_hydrogen_header(arg[1][-8:-4], "Hydrogen bonds main-main", [3.5, 4])
    pos_prev = []
    mmhb = parsing(arg[1], ["N", "O", "OXT"], [])
    for i, elem1 in enumerate(mmhb):
        for elem2 in mmhb[i + 1:]:
            # check the donor-acceptor pairs
            if elem1.name == "N" and elem2.name in ["O", "OXT"]:
                donor = elem1
                acceptor = elem2
            elif elem2.name == "N" and elem1.name in ["O", "OXT"]:
                donor = elem2
                acceptor = elem1
            else:
                continue
            dist = calc_range(donor, acceptor)
            if check_criteria(arg[2], donor, acceptor, pos_prev, dist, 3.5):
                if abs(donor.position - acceptor.position) >= 2 and donor.residue != "PRO":
                    print_hydrogen_res(arg[1][-8:-4], donor.position, donor.residue, donor.chain,
                                       donor.name, acceptor.position, acceptor.residue,
                                       acceptor.chain, acceptor.name, dist)
                    pos_prev.append([donor.position, acceptor.position,
                                     donor.chain, acceptor.chain])
    print("\n")


def ms_hbond(arg):
    """
    function calculating the main chain-side chain hydrogen bonds

    :param arg: the list of arguments of the command line
    """
    print_hydrogen_header(arg[1][-8:-4], "Hydrogen bonds main-side", [3.5, 4])
    main = ["N", "O", "OXT"]
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
    mshb = parsing(arg[1], main + side_all, [])
    for i, elem1 in enumerate(mshb):
        for elem2 in mshb[i + 1:]:
            # check the donor-acceptor pairs
            def_range = 3.5
            if elem1.name in main and elem2.name in side_all:
                if elem1.name == "N" and \
                        (elem2.residue in side_acc.keys() and
                         elem2.name in side_acc[elem2.residue]):
                    donor = elem1
                    acceptor = elem2
                elif (elem1.name == "O" or elem1.name == "OXT") and \
                        (elem2.residue in side_don.keys() and
                         elem2.name in side_don[elem2.residue]):
                    donor = elem2
                    acceptor = elem1
            elif elem1.name in side_all and elem2.name in main:
                if elem2.name == "N" and \
                        (elem1.residue in side_acc.keys() and
                         elem1.name in side_acc[elem1.residue]):
                    donor = elem2
                    acceptor = elem1
                elif (elem2.name == "O" or elem2.name == "OXT") and \
                        (elem1.residue in side_don.keys() and
                         elem1.name in side_don[elem1.residue]):
                    donor = elem1
                    acceptor = elem2
            else:
                continue
            if donor is not None and acceptor is not None and (donor, acceptor) != prev_do_ac:
                if donor.residue == "PRO" and donor.name == "N":
                    continue
                if donor.name == "SG" or acceptor.name == "SD":
                    def_range = 4
                dist = calc_range(donor, acceptor)
                if dist <= def_range and \
                        ([donor.position, acceptor.position, donor.chain, acceptor.chain]
                         not in pos_prev or donor.name != acceptor.name) and \
                        donor.check_non_id_hbond_side(acceptor) and \
                        abs(donor.position - acceptor.position) >= 2:
                    if arg[2] == "--intra" and elem1.chain == elem2.chain:
                        print_hydrogen_res(arg[1][-8:-4], donor.position, donor.residue,
                                           donor.chain, donor.name, acceptor.position,
                                           acceptor.residue, acceptor.chain, acceptor.name, dist)
                        pos_prev.append([donor.position, acceptor.position, donor.chain,
                                         acceptor.chain])
                    elif arg[2] == "--inter" and elem1.chain != elem2.chain:
                        print_hydrogen_res(arg[1][-8:-4], donor.position, donor.residue,
                                           donor.chain, donor.name, acceptor.position,
                                           acceptor.residue, acceptor.chain, acceptor.name, dist)
                        pos_prev.append([donor.position, acceptor.position, donor.chain,
                                         acceptor.chain])
                    else:
                        continue
                prev_do_ac = (donor, acceptor)
    print("\n")


def ss_hbond(arg):
    """
    function calculating the main side-side chain hydrogen bonds

    :param arg: the list of arguments of the command line
    """
    print_hydrogen_header(arg[1][-8:-4], "Hydrogen bonds side-side", [3.5, 4])
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
    sshb = parsing(arg[1], side_all, [])
    for elem1 in sshb:
        for elem2 in sshb:
            # check the donor-acceptor pairs
            def_range = 3.5
            if (elem1.residue in poss_don.keys() and elem1.name in poss_don[elem1.residue]) and \
                    (elem2.residue in poss_acc.keys() and elem2.name in poss_acc[elem2.residue]):
                donor = elem1
                acceptor = elem2
            elif (elem1.residue in poss_acc.keys() and elem1.name in poss_acc[elem1.residue]) and \
                    (elem2.residue in poss_don.keys() and elem2.name in poss_don[elem2.residue]):
                donor = elem2
                acceptor = elem1
            else:
                continue
            if donor.name == "SG" or acceptor.name == "SD":
                def_range = 4
            dist = calc_range(donor, acceptor)
            if dist <= def_range and \
                    ([donor.position, acceptor.position, donor.chain, acceptor.chain]
                     not in pos_prev) and \
                    donor.check_non_id_hbond_side(acceptor) and \
                    abs(donor.position - acceptor.position) >= 2:
                if arg[2] == "--intra" and elem1.chain == elem2.chain:
                    print_hydrogen_res(arg[1][-8:-4], donor.position, donor.residue, donor.chain,
                                       donor.name, acceptor.position, acceptor.residue,
                                       acceptor.chain, acceptor.name, dist)
                    pos_prev.append([donor.position, acceptor.position,
                                     donor.chain, acceptor.chain])
                elif arg[2] == "--inter" and elem1.chain != elem2.chain:
                    print_hydrogen_res(arg[1][-8:-4], donor.position, donor.residue, donor.chain,
                                       donor.name, acceptor.position, acceptor.residue,
                                       acceptor.chain, acceptor.name, dist)
                    pos_prev.append([donor.position, acceptor.position,
                                     donor.chain, acceptor.chain])
                else:
                    continue
    print("\n")


def print_header(pdb_name, title, range):
    """
    function printing and writing in the res file the head of the result's table
    """
    with open("../results/"+pdb_name+"_res.txt", "a") as fout:
        # write in the result file
        if title == "Aromatic-aromatic interactions":
            print("\n\n{:>45} (between {}-{}A)".format(title, range[0], range[1]))
            fout.write("\n\n\n{:>45} (between {}-{}A)".format(title, range[0], range[1]))
        else:
            print("\n\n{:>45} (<= {}A)".format(title, range))
            fout.write("\n\n\n{:>45} (<= {}A)".format(title, range))
        print("{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}".format("Position",
                                                                  "Residue", "Chain",
                                                                  "Position", "Residue",
                                                                  "Chain", "Distance"))
        fout.write("\n{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}".format("Position",
                                                                         "Residue", "Chain",
                                                                         "Position", "Residue",
                                                                         "Chain", "Distance"))


def print_pos_res_ch_dis(pdb_name, pos1, res1, chain1, pos2, res2, chain2, dist):
    """
    function printing and writing in the res file the result's table
    """
    with open("../results/" + pdb_name + "_res.txt", "a") as fout:
        # write in the result file
        ind = 0
        print("{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10.2f}"
              .format(pos1, res1, chain1, pos2, res2, chain2, dist))
        if ind == 0:
            # if it's the first line of results
            fout.write("\n{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10.2f}"
                       .format(pos1, res1, chain1, pos2, res2, chain2, dist))
            ind += 1
        else:
            fout.write("{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10.2f}"
                       .format(pos1, res1, chain1, pos2, res2, chain2, dist))


def print_hydrogen_header(pdb_name, title, range):
    """
    function printing and writing in the res file the head of the result's table
    """
    with open("../results/"+pdb_name+"_res.txt", "a") as fout:
        # write in the result file
        print("\n\n{:>45} (<= {}/{}A)".format(title, range[0], range[1]))
        fout.write("\n\n\n{:>45} (<= {}/{}A)".format(title, range[0], range[1]))
        print("{:^10}{:^10}{:^10}{:^10}{:^10}"
              "{:^10}{:^10}{:^10}{:^10}".format("Position", "Residue", "Chain",
                                                "Atom", "Position", "Residue", "Chain",
                                                "Atom", "Distance"))
        fout.write("\n{:^10}{:^10}{:^10}{:^10}{:^10}"
                   "{:^10}{:^10}{:^10}{:^10}".format("Position", "Residue", "Chain",
                                                     "Atom", "Position", "Residue",
                                                     "Chain", "Atom", "Distance"))


def print_hydrogen_res(pdb_name, posd, resd, chaind, named, posa, resa, chaina, namea, dist):
    """
    function printing and writing in the res file the the result's table
    """
    with open("../results/" + pdb_name + "_res.txt", "a") as fout:
        # write in the result file
        ind = 0
        print("{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10.2f}"
              .format(posd, resd, chaind, named, posa, resa, chaina, namea, dist))
        if ind == 0:
            # if it's the first line of results
            fout.write("\n{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10.2f}"
                       .format(posd, resd, chaind, named, posa, resa, chaina, namea, dist))
            ind += 1
        else:
            fout.write("{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10}{:^10.2f}"
                       .format(posd, resd, chaind, named, posa, resa, chaina, namea, dist))
