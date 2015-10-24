# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 16:31:41 2015

@author: jgb
"""

import os, importlib, shutil
import numpy as np
import numpy.random as rand
from numpy import pi, sin, cos, sqrt

def jackknife_var_on_mean(array, block_size=1):
    n = np.size(array)
    mean = np.mean(array)
    var = np.var(array)
    var_mean = var/(n)
    
    jack_ests = np.array([(1/(n-1))*np.sum(np.concatenate((array[:i],array[i+1:]))) for i in range(n)])
    jack_var_mean = ((n-1)/float(n))*np.sum((jack_ests-mean)**2)
    print("Std var mean: {:.3g}, Jackknife var mean: {:.3g}".format(var_mean,jack_var_mean))
    
    return jack_var_mean
    
def jackknife_err_on_var(array, block_size=1):
    n = np.size(array)
    var = np.var(array)
    
    jack_pseudos = np.array([n*var - (n-1)*np.var(np.concatenate((array[:i],array[i+1:]))) for i in range(n)])
    jack_pseudo_mean = (1/float(n))*np.sum(jack_pseudos)  
    jack_var_var = (1/float(n-1))*np.sum((jack_pseudos-jack_pseudo_mean)**2)
    print("Std err var: {:.3g}, Jackknife err var: {:.3g}".format(sqrt(var)/sqrt(2*(n-1)),sqrt(jack_var_var)/2.0/sqrt(n)))
    
    return jack_pseudo_mean,jack_var_var

#tests    
if __name__=='__main__':
    data = rand.normal(loc=0.0, scale=1.0, size=500)
    jack_var_mean = jackknife_var_on_mean(data)  
    frac_diff = jack_var_mean/(np.var(data)/np.size(data))
    print(frac_diff)
    if frac_diff>0.99 and frac_diff<1.01:
        print("Jacknife variance test with Gaussian: OK")
    else:
        print("Jacknife variance test with Gaussian: FAIL")
        
    jack_pseudo_mean,jack_var_var = jackknife_err_on_var(data)
    print(jack_pseudo_mean)