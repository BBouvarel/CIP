#!/usr/bin/env python

class Atom:
    """
    class representing the atoms and some of their characteristics contained in a pdb file

    var name: atom name
        element: element of the atom
        position: position of the residue that contains the atom
        residue: the residue containing the atom
        chain: the chains containing the atom
        x: coordinate x of the atom
        y: coordinate y of the atom
        z: coordinate z of the atom
    """

    def __init__(self, name, element, position, residue, chain, x, y, z):
        self.name = name
        self.element = element
        self.position = int(position)
        self.residue = residue
        self.chain = chain
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)