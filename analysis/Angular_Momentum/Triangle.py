"""

Tests the triangle inequalities in Messiah, A. "Quantum Mechanics" v. 2,
pg. 1062 (North-Holland, Amsterdam) 1962.  Returns True if the inequalities
are satisfied, false if not.  Also tests if the triad sums to an integer

Written:  KAE University at Albany Physics Department 26 Oct 08

"""
from numpy import *

def Triangle(x,y,z):
    if ((abs(x-y) <= z) and (z <= x+y) and (floor(x+y+z) == x+y+z)):
        test = True
    else:
        test = False

    return test
