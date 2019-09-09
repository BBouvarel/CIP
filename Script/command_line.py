#!/usr/bin/env python

import sys

def check_command(argv):
    """
    function to detect errors in the command line

    param argv: list of the elements of the command line
    return: a value depending of the error
    """
    if len(argv) < 2:
        return 1
    if argv[1][-4:] != ".pdb":
        return 2
    if len(argv) == 2:
        return 3

def set_val_default(arg_input, arg_name, def_range):
    """
    function to assign a threshold distance value used by some calculation of interactions

    param arg_input: argument calling a specific interaction calculation write by the user
                     (can be associate to a number)
          arg_name: name of the argument (without number)
          def_range: default value used by the interaction calculation
    return: a threshold distance or -1 if there is an error in the number send by the user
    """
    range = def_range
    if len(arg_input) > len(arg_name):
        # while True:
        try:
            float(arg_input[len(arg_name):])
            range = float(arg_input[len(arg_name):])
            # break
        except ValueError:
            print("The value of the argument", arg_name, "is invalid")
            range = -1
    return range

def check_help(argv):
    """
    function to check the presence of the argument help and to run it

    param argv: list of the elements of the command line
    return: a boolean to exit the program
    """
    if len(argv) == 2 and argv[0][-9:] == "script.py" and argv[1] == "--help":
        print("\nHelp:\n\n"
              "Command:\n"
              "python script.py file.pdb --arg1 --arg2 --argN\n\n"
              "Options:\n"
              "--hphb  ->  Run the calculation of hydrophobic interactions.\n"
              "            --hphbVALUE Give a specific distance value.\n"
              "            (ex: --hphb1.0 , default value: 5A).\n\n"
              "--inic  ->  Run the calculation of ionic interactions.\n"
              "            --inicVALUE Give a specific distance value.\n"
              "            (ex: --inic1.0 , default value: 6A).\n\n"
              "--arar  ->  Run the calculation of aromatic-aromatic interactions.\n"
              "            --ararVALUE Give a specific distance value.\n"
              "            (ex: --arar1.0 , default value: 4.5A).\n\n"
              "--arsu  ->  Run the calculation of aromatic-sulphur interactions.\n"
              "            --arsuVALUE Give a specific distance value.\n"
              "            (ex: --arsu1.0 , default value: 5.3A).\n\n"
              "--capi  ->  Run the calculation of cation-pi interactions.\n"
              "            --capiVALUE Give a specific distance value.\n"
              "            (ex: --capi1.0 , default value: 6A).\n\n"
              "--disu  ->  Run the calculation of disulphide bridges.\n\n"
              "--mmhb  ->  Run the calculation of main chain-main chain hydrogen bond.\n\n"
              "--mshb  ->  Run the calculation of main chain-side chain hydrogen bond.\n\n"
              "--sshb  ->  Run the calculation of side chain-side chain hydrogen bond.\n\n"
              )
        return True