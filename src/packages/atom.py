#!/usr/bin/env python

"""
package containing the atom class in which are stored the atoms of
the PDB file that meet the criteria for a specific interaction
"""

__author__ = "Bertrand BOUVAREL"
__date__ = "2019/09"


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
    def __init__(self, name, position, residue, chain, coor_x, coor_y, coor_z):
        self.name = name
        self.position = int(position)
        self.residue = residue
        self.chain = chain
        self.coor_x = float(coor_x)
        self.coor_y = float(coor_y)
        self.coor_z = float(coor_z)

    def check_non_id(self, ato2):
        """
        function checking if the atom in the argument is identical to the this atom

        :param ato2: object of the class atom compared
        :return: boolean to know if the two atoms are identical
        """
        if self.position != ato2.position or self.chain != ato2.chain:
            return True
        return False

    def check_non_id_hbond_side(self, ato2):
        """
        function checking if the atom in the argument is identical to the this atom
        version for m-s and s-s hbond

        :param ato2: object of the class atom compared
        :return: boolean to know if the two atoms are identical
        """
        if self.position != ato2.position or self.chain != ato2.chain or self.name != ato2.name:
            return True
        return False
