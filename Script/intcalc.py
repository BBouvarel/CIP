#!/usr/bin/env python

import math

def calc_range3D(elem1, elem2):
    """
    function allowing the calculation of the distance in three dimensions of two atoms

    param elem1, elem2: two objects of the class atom
    return: distance between two atoms
    """
    range = math.sqrt((elem2.x-elem1.x)**2+(elem2.y-elem1.y)**2+(elem2.z-elem1.z)**2)
    return range
