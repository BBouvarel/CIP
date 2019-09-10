#!/usr/bin/env python

import math
import atom


def calc_range3D(elem1, elem2):
    """
    function allowing the calculation of the distance in three dimensions of two atoms

    param elem1, elem2: two objects of the class atom
    return: distance between two atoms
    """
    dist = math.sqrt((elem2.x-elem1.x)**2+(elem2.y-elem1.y)**2+(elem2.z-elem1.z)**2)
    return dist

def add_to_list(line, liste):
    liste.append(atom.Atom(line[12:15].strip(), line[76:78].strip(),
                              line[22:26].strip(), line[17:20].strip(),
                              line[21:22].strip(), line[30:38].strip(),
                              line[38:46].strip(), line[46:54].strip()))
    return liste

def parsing(pdb, name, res, symbol):
    """
    function allowing the parsing of the pdb file following certain conditions

    param pdb: the pdb file given by the user
          name: atom name
          res: list of residues containing the element
          symbol: list of elements symbols
    return: list of objects of class atom that respect the conditions
    """
    elements = []
    #if len(name) == 0 : continue
    #if len(res) == 0: res = [""]
    #if len(symbol) == 0: symbol = ["C", "O", "N", "S"]

    # juste mettre if des conditions l'un apres l autre avec and len(liste) != 0

    # if len res symbl ou autre =0 alors on initialise ces liste avec toutes les possibilités
    with open(pdb, "r") as finput:
        for line in finput:
            if line[0:6].strip() == "ATOM":
                if line[17:20].strip() in res and len(res) != 0:
                    if line[76:78].strip() in symbol and len(symbol) != 0:
                        elements = add_to_list(line, elements)
                    elif line[12:15].strip() in name and len(name) != 0:
                        # if no restriction on the element's symbol
                        elements = add_to_list(line, elements)
                            #elements.append(atom.Atom(line[12:16].strip(), line[76:78].strip(),
                                                      #line[22:26].strip(), line[17:20].strip(),
                                                      #line[21:22].strip(), line[30:38].strip(),
                                                      #line[38:46].strip(), line[46:54].strip()))
    return elements