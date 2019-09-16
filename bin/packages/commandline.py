#!/usr/bin/env python

"""
package for checking errors in the command line passed to cip.py,
and assigning threshold distance values if necessary
"""

__author__ = "Bertrand BOUVAREL"
__date__ = "2019/09"


def check_command(argv):
    """
    function to detect errors in the command line

    :param argv: list of the elements of the command line
    :return: an int depending of the error
    """
    if len(argv) < 2:
        return 1
    if argv[1][-4:] != ".pdb":
        # use of a pdb file
        return 2
    if len(argv) == 2:
        # no arguments for the interactions calculation
        return 3
    if argv[2] not in ["--intra", "--inter"]:
        return 4


def set_val_default(arg_input, arg_name, def_range):
    """
    function to assign a threshold distance value used by some calculation of interactions

    :param arg_input: argument calling a specific interaction calculation write by the user
                      (can be associate to a number)
    :param arg_name: name of the argument (without number)
    :param def_range: default value used by the interaction calculation
    :return: a threshold distance or -1 if there is an error in the number send by the user
    """
    range = def_range
    if len(arg_input) > len(arg_name):
        if arg_name == "--arar":
            # two distance threshold values to make an interval
            try:
                # check if the two values are numbers and not a letters or a symbols
                val1 = arg_input[len(arg_name):].split("/")[0]
                val2 = arg_input[len(arg_name):].split("/")[1]
                if val1 >= val2:
                    range = -1
                    print("The values with the argument", arg_name,
                          "are invalid, please leave a positive gap between the two values")
                else:
                    range = (float(val1), float(val2))
            except ValueError:
                print("The values with the argument", arg_name, "are invalid")
                range = -1
        else:
            # only one threshold value needed
            try:
                # check if the value with the arg is a number and not a letter or a symbol
                range = float(arg_input[len(arg_name):])
            except ValueError:
                print("The value with the argument", arg_name, "is invalid")
                range = -1
    return range


def check_help(argv):
    """
    function to check the presence of the argument help and to run it

    :param argv: list of the elements of the command line
    :return: a boolean to exit the program
    """
    if len(argv) == 2 and argv[0][-6:] == "cip.py" and argv[1] == "--help":
        print("\nHelp:\n\n"
              "Command:\n"
              "python3 cip.py ../data/file.pdb --intra --arg1 --arg2 --argN\n\n"
              "Options:\n"
              "--help  ->  Show the help message\n\n"
              "--intra / --inter ->  calculates intra / inter protein interactions\n\n"
              "--hphb  ->  Run the calculation of hydrophobic interactions\n"
              "            --hphbVALUE Give a specific distance value\n"
              "            (ex: --hphb1.0 , default value: 5A)\n\n"
              "--inic  ->  Run the calculation of ionic interactions\n"
              "            --inicVALUE Give a specific distance value\n"
              "            (ex: --inic1.0 , default value: 6A)\n\n"
              "--arar  ->  Run the calculation of aromatic-aromatic interactions\n"
              "            --ararVALUE1/VALUE2 Give a specific distance interval\n"
              "            (ex: --arar1.0/5.0 , default value: 4.5A to 7A)\n\n"
              "--arsu  ->  Run the calculation of aromatic-sulphur interactions\n"
              "            --arsuVALUE Give a specific distance value\n"
              "            (ex: --arsu1.0 , default value: 5.3A)\n\n"
              "--capi  ->  Run the calculation of cation-pi interactions\n"
              "            --capiVALUE Give a specific distance value\n"
              "            (ex: --capi1.0 , default value: 6A)\n\n"
              "--disu  ->  Run the calculation of disulphide bridges\n\n"
              "--mmhb  ->  Run the calculation of main chain-main chain hydrogen bond\n\n"
              "--mshb  ->  Run the calculation of main chain-side chain hydrogen bond\n\n"
              "--sshb  ->  Run the calculation of side chain-side chain hydrogen bond\n\n"
              "Examples of command:\n"
              "> python3 cip.py ../data/1atn.pdb --intra --inic --disu --arar --arsu --capi "
              "--hphb --mmhb --mshb --sshb\n\n"
              "> python3 cip.py ../data/1atn.pdbv --inter --inic --disu --arar --arsu --capi "
              "--hphb --mmhb --mshb --sshb\n\n"
              "> python3 cip.py ../data/1atn.pdb --intra --inic6.5 --disu --arar4/6.5 --arsu5.5 "
              "--capi7 --hphb4.5 --mmhb --mshb --sshb\n\n")
        return True
