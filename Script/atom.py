#!/usr/bin/env python

class Atom:
    """
    class representing the atoms and some of their characteristics contained in a pdb file

    :param name: atom name
    :param element: element of the atom
    :param position: position of the residue that contains the atom
    :param residue: the residue containing the atom
    :param chain: the chains containing the atom
    :param x: coordinate x of the atom
    :param y: coordinate y of the atom
    :param z: coordinate z of the atom
    """

    def __init__(self, name, position, residue, chain, x, y, z):
        self.name = name
        self.position = int(position)
        self.residue = residue
        self.chain = chain
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def get_name(self):
        """
        function returning the atom name
        """
        return self.name

    def get_position(self):
        """
        function returning the position of the atom
        """
        return self.position

    def get_residue(self):
        """
        function returning the residue which contains the atom
        """
        return self.residue

    def get_chain(self):
        """
        function returning the chain which contains the atom
        """
        return self.chain

    def get_x(self):
        """
        function returning the coordinate x of the atom
        """
        return self.x

    def get_y(self):
        """
        function returning the coordinate y of the atom
        """
        return self.y

    def get_z(self):
        """
        function returning the coordinate z of the atom
        """
        return self.z
