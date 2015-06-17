"""
Python implementation of Lai's 6j program originally developed for Maple.
As nearly as I can tell, there are no checks of the triangle conditions, so
it's up to the user to verify that all arguments lead to an allowed 3j symbol.

The output has been tested against Anthony Stone's Wigner Coefficient Calculator

http://www-stone.ch.cam.ac.uk/wigner.html

for a few cases.  More extensive testing would be needed to ensure that this
module was road worthy.

Literature: Lai, S.-T. "Computation of Algebraic Formulas for Wigner 3-j, 6-j,
and 9-j Symbols by Maple" International Journal of Quantum Chemistry,
52:593--607 (1994)

Arguments are called in the order

           { e a f }
           { b d c }

according to the notation in Lai's paper

Written:  KAE Department of Physics, University at Albany 26 Oct 08
Uses:     numpy, scipy
Calls:    Delta (sqrt of ratios of factorials), Triangle (tests triangle
          conditions)
Returns:  Wigner 6J coefficient in s1, flag True if all constraints satisfied
Future:   Explicitly check triangle conditions
Modified: 27 Oct 08 implement future

"""

from numpy import *
from scipy import *

from Delta import *
from Triangle import *

def W6J(e,a,f,b,d,c):
    if (Triangle(e,a,f) and Triangle(e,d,c) and Triangle(b,a,c) and \
        Triangle(b,d,f)):
        flag = True
    else:
        flag = False

    if flag:
        aa = factorial(a+b+c+1)*factorial(b+d+f+1)/\
             factorial(a+b-c)/factorial(c-d+e)/factorial(c+d-e)/\
             factorial(-e+a+f)/factorial(e-a+f)/factorial(b+d-f)
        v  = min(2*b,-a+b+c,b-d+f)
        s0 = 0.
        for n in arange(0,v+1):
            s0 = s0 +(-1)**n*factorial(2*b-n)/factorial(n)* \
                 factorial(b+c-e+f-n)/factorial(-a+b+c-n)* \
                 factorial(b+c+e+f+1-n)/factorial(b-d+f-n)/\
                 factorial(a+b+c+1-n)/factorial(b+d+f+1-n)

        s1 = ((-1)**(b+c+e+f))*Delta(a,b,c)*Delta(c,d,e)*Delta(a,e,f)* \
             Delta(b,d,f)*aa*s0

    else:
        s1 = 0.

    return s1,flag
