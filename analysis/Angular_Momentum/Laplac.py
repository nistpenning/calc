from numpy import *
from scipy import *

def Laplac(a,b,c):
    l = sqrt(factorial(a-b+c)*factorial(a+b-c)*factorial(a+b+c+1)/\
             factorial(-a+b+c))
    return l
