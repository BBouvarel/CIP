#!/usr/bin/env python

import sys

def check_command(argv):
    if len(argv) < 2:
        return 1
    if argv[1][-4:] != ".pdb":
        return 2
    if len(argv) == 2:
        return 3

def set_val_default(arg_input, arg_name, range):
    def_range = range
    if len(arg_input) > len(arg_name):
        # while True:
        try:
            float(arg_input[len(arg_name):])
            def_range = float(arg_input[len(arg_name):])
            # break
        except ValueError:
            print("La valeur de l'argument", arg_name, "est invalide")
            def_range = -1
    return def_range

def check_help(argv):
    if len(argv) == 2 and argv[0][-9:] == "script.py" and argv[1] == "--help":
        print("HELP\n\n"
              "Commandes:\n"
              "--hphb  ->  lance le calcul d'interactions hydrophobes.\n"
              "            --hphbNOMBRE donne une valeur de distance\n"
              "            (ex: --hphb2 , valeur par defaut: 5A).\n\n"
              "--inic  ->  lance le calcul d'interactions ioniques.\n"
              "            --inicNOMBRE donne une valeur de distance\n"
              "            (ex: --inic2 , valeur par defaut: 6A).\n\n"
              )
        return True