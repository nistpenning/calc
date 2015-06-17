# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 12:35:20 2014

@author: bsawyer
"""
from __future__ import division

from numpy import *
from scipy import *
from math import *

#from Delta import *
#from Triangle import *
#from Laplac import *

def W3J(J1,J,J2,M1,M,M2):

    # set flags to determine allowed 3j symbols based on calling arguments
    
    if ((abs(J1-J) <= J2) and (J2 <= J1+J)):
        tri = True
    else:
        tri = False

    if ((M1 in arange(-J1, J1+1)) and (M in arange(-J, J+1)) and \
        (M2 in arange(-J2, J2+1))):
        mem = True
    else:
        mem = False

    if (floor(J1 + J + J2) == J1 + J + J2):
        perim = True
    else:
        perim = False

    if (M1 + M + M2 == 0):
        mag = True
    else:
        mag = False

    # now compute allowed 3j symbol, return 0 if not allowed

    flag = tri and mem and perim and mag
    
    if flag:
        delta = sqrt(float(factorial(J1+J2-J)*float(factorial(J1-J2+J))* \
                     float(factorial(-J1+J2+J)))/float(factorial(J1+J+J2+1)))

        a = sqrt(float(factorial(J2+M2)*float(factorial(J2-M2)))/float(factorial(J+M))/\
                 float(factorial(J-M))/float(factorial(J1-M-M2))/float(factorial(J1+M+M2)))
                           
        s0 = 0.0

        u  = -J1+J2+J
                
        z1 = max(-J1-M-M2,u-J2-M2,0)
        z2 = min(J+J2-M-M2,J2-M2,u)
        
        for z in arange(z1,z2+1):
            if ((0 <= J+J2-M-M2-z) and (0 <= J2-M2-z) and \
                (0 <= u-z) and (0 <= J2+M2-u+z)):
                s0 = float(s0) + (-1.0)**(2*J-J1-M1+z) * \
                     float(factorial(J+J2-M-M2-z))*float(factorial(J1+M+M2+z))/\
                     float(factorial(z))/float(factorial(J2-M2-z))/float(factorial(u-z))/\
                     float(factorial(J2+M2-u+z))

        s = a*delta*s0

    else:
        s = 0

    return s


def W6J(e,a,f,b,d,c):
    if (Triangle(e,a,f) and Triangle(e,d,c) and Triangle(b,a,c) and \
        Triangle(b,d,f)):
        flag = True
    else:
        flag = False

    if flag:
        aa = float(factorial(a+b+c+1.0))*float(factorial(b+d+f+1.0))*1.0/\
             factorial(a+b-c)/factorial(c-d+e)*1.0/factorial(c+d-e)*1.0/\
             factorial(-1.0*e+a+f)/factorial(e-a+f)/factorial(b+d-f)
        v  = min(2*b,-a+b+c,b-d+f)
        s0 = 0.0
        for n in arange(0,v+1):
            s0 = s0 +(-1.0)**n*float(factorial(2*b-n))/factorial(n)*1.0* \
                 factorial(b+c-e+f-n)/factorial(-a+b+c-n)* \
                 factorial(b+c+e+f+1.0-n)/factorial(b-d+f-n)/\
                 factorial(a+b+c+1.0-n)/factorial(b+d+f+1.0-n)

        s1 = ((-1.0)**(b+c+e+f))*Delta(a,b,c)*Delta(c,d,e)*Delta(a,e,f)*1.0* \
             Delta(b,d,f)*aa*s0

    else:
        s1 = 0.0

    return s1


def W9J(j11,j21,j31,j12,j22,j32,j13,j23,j33):
    
    if (floor(j13+j23-j33) == j13+j23-j33):
        s0 = (-1.0)**(j13+j23-j33)*Laplac(j21,j11,j31)*Laplac(j12,j22,j32)/\
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
                                      factorial(2.0*j23-x)*
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
    
    return s2

def Laplac(a,b,c):
    l = sqrt(factorial(a-b+c)*factorial(a+b-c)*factorial(a+b+c+1.0)/\
             factorial(-a+b+c))
    return l

def Triangle(x,y,z):
    if ((abs(x-y) <= z) and (z <= x+y) and (floor(x+y+z) == x+y+z)):
        test = True
    else:
        test = False

    return test

def Delta(a,b,c):
    d = sqrt(factorial(a+b-c)*factorial(a-b+c)/factorial(a+b+c+1.0)*\
             factorial(-a+b+c))
    return d

