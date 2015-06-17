from numpy import *
from scipy import *

def Delta(a,b,c):
    d = sqrt(factorial(a+b-c)*factorial(a-b+c)/factorial(a+b+c+1)*\
             factorial(-a+b+c))
    return d
