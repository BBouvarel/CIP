#!/usr/bin/env python

"""
    program used to launch the intra-protein interaction calculations by
    checking the conformity of the command line (with the commandline.py
    package), then sending to intcalc.py the desired interactions
"""

__author__ = "Bertrand Bouvarel"
__date__ = "2019/09"


import sys
import packages.commandline as cl
import packages.intcalc as ic


if __name__ == "__main__":
    # interaction using distance threshold and the default values
    INTERAC_WITH_RANGE = ["--hphb", "--inic", "--arar", "--arsu", "--capi"]
    INTERAC_RANGE = [5, 6, [4.5, 7], 5.3, 6]

    if cl.check_help(sys.argv):
        # check if the help is called
        sys.exit()
    # check for potential errors in the command line
    if cl.check_command(sys.argv) == 1:
        sys.exit("A pdb file or arguments are missing, type --help to access the program help")
    elif cl.check_command(sys.argv) == 2:
        sys.exit("A pdb file is needed, type --help to access the program help")
    elif cl.check_command(sys.argv) == 3:
        sys.exit("No interactions calculation are given, type --help to access the program help")

    with open("../results/"+sys.argv[1][-8:-4]+"_res.txt", "w") as fout:
        fout.write("Results of the intra-protein interaction calculation:")
        # create an empty result file

    for arg in range(2, len(sys.argv)):
        # run all the calculation arguments of the command line
        if sys.argv[arg][0:6] in INTERAC_WITH_RANGE:
            # check if the calculation of the interaction need a distance threshold
            def_range = cl.set_val_default(sys.argv[arg], sys.argv[arg][0:6],
                                           INTERAC_RANGE[INTERAC_WITH_RANGE.index(
                                               sys.argv[arg][0:6])])
            if def_range == -1:
                # error in the command line threshold value, the program skip the argument
                continue

        if sys.argv[arg][0:6] == "--hphb":
            # calculation of hydrophobic interaction
            ic.hydrophobic(sys.argv, def_range)

        elif sys.argv[arg][0:6] == "--inic":
            # calculation of ionic interaction
            ic.ionic(sys.argv, def_range)

        elif sys.argv[arg][0:6] == "--arar":
            # calculation of aromatic-aromatic interaction
            ic.aro_aro(sys.argv, def_range)

        elif sys.argv[arg][0:6] == "--arsu":
            # calculation of aromatic-sulphur interaction
            ic.aro_sul(sys.argv, def_range)

        elif sys.argv[arg][0:6] == "--capi":
            # calculation of Cation-pi interaction
            ic.cation_pi(sys.argv, def_range)

        elif sys.argv[arg][0:6] == "--disu":
            # calculation of disulphide bridges
            ic.disulphide(sys.argv)

        elif sys.argv[arg][0:6] == "--mmhb":
            # calculation of hydrogen bonds main-main
            ic.mm_hbond(sys.argv)

        elif sys.argv[arg][0:6] == "--mshb":
            # calculation of hydrogen bonds main-side
            ic.ms_hbond(sys.argv)

        elif sys.argv[arg][0:6] == "--sshb":
            # calculation of hydrogen bonds side-side
            ic.ss_hbond(sys.argv)

        else:
            sys.exit("One of the arguments is not recognized by the program, "
                     "type --help to access the program help")
