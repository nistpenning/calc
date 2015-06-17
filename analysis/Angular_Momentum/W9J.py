"""
Python implementation of Lai's 9j program originally developed for Maple.
As nearly as I can tell, there are no checks of the triangle conditions, so
it's up to the user to verify that all arguments lead to an allowed 9j symbol.

The output has been tested against Anthony Stone's Wigner Coefficient Calculator

http://www-stone.ch.cam.ac.uk/wigner.html

for a few cases.  More extensive testing would be needed to ensure that this
module was road worthy.

The arguments are ordered as follows:

            (j11 j21 j31)
            (j12 j22 j32)
            (j13 j23 j33)

in Lai's notation.

Literature: Lai, S.-T. "Computation of Algebraic Formulas for Wigner 3-j, 6-j,
and 9-j Symbols by Maple" International Journal of Quantum Chemistry,
52:593--607 (1994)

Written:  KAE Department of Physics, University at Albany 26 Oct 08
Uses:     numpy, scipy
Calls:    Laplac (sqrt of ratios of factorials)
Returns:  Wigner 9J coefficient in s2, flag is True if calling arguments OK.
Future:   Explicitly check factorial and phase arguments
Modified: Enforced stricter checks for non-negative factorial arguments (not
          in Lai's codes) 27 Oct 08

"""

from numpy import *
from scipy import *

from Laplac import *

def W9J(j11,j21,j31,j12,j22,j32,j13,j23,j33):
    
    if (floor(j13+j23-j33) == j13+j23-j33):
        s0 = (-1)**(j13+j23-j33)*Laplac(j21,j11,j31)*Laplac(j12,j22,j32)/\
             Laplac(j21,j22,j23)/Laplac(j12,j11,j13)*Laplac(j33,j31,j32)/\
             Laplac(j33,j13,j23)
        sflag = True
    else:
        s0 = 0
        sflag = False

    u1 = min(2*j23,j22-j21+j23,j13+j23-j33)
    if (u1 >= 0):
        xflag = True
    else:
        xflag = False

    v1 = min(j12+j22-j32,j31-j32+j33)
    if (v1 >= 0):
        yflag = True
    else:
        yflag = False

    w1 = min(j11-j12+j13,j11+j21+j31+1,j11+j21-j31,2*j11)
    if (w1 >= 0):
        zflag = True
    else:
        zflag = False

    
    s1 = 0.

    if (sflag and xflag and yflag and zflag):
        for x in arange(0, u1+1):
            if (((j21+j22-j23+x) >= 0) and ((j13-j23+j33+x) >= 0)):
                for y in arange(0, v1+1):
                    if (((j22-j12+j32+y) >= 0) and ((j31+j32-j33+y) >= 0)
                        and ((j21-j12+j32-j23+x+y) >= 0)):
                        for z in arange(0, w1+1):
                            if (((j12-j11+j13) >= 0) and
                                ((j12-j11-j23+j33+x+z) >= 0) and
                                ((j11+j21-j32+j33-y-z) >= 0)):
                                s1 = (s1+(-1.)**(x+y+z)*
                                      factorial(2*j23-x)*
                                      factorial(j21+j22-j23+x)/
                                      factorial(x)/
                                      factorial(j22-j21+j23-x)/
                                      factorial(j13+j23-j33-x)/
                                      factorial(j21-j12+j32-j23+x+y)*
                                      factorial(j13-j23+j33+x)*
                                      factorial(j22-j12+j32+y)/
                                      factorial(-j11+j12+j33-j23+x+z)/
                                      factorial(y)/
                                      factorial(j12+j22-j32-y)*
                                      factorial(j31+j32-j33+y)*
                                      factorial(j11+j21-j32+j33-y-z)/
                                      factorial(j31-j32+j33-y)/
                                      factorial(2*j32+1+y)/
                                      factorial(z)/
                                      factorial(j11+j21-j31-z)*
                                      factorial(2*j11-z)*
                                      factorial(j12-j11+j13+z)/
                                      factorial(j11-j12+j13-z)/
                                      factorial(j11+j21+j31+1-z))
                            else:
                                continue
                    else:
                        continue
            else:
                continue
        flag = True
        s2 = s0*s1
    else:
        flag = False
        s2 = 0
    
    return s2,flag


