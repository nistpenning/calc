# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 16:31:41 2015

@author: jgb
"""

import os, importlib, shutil
import numpy as np
import numpy.random as rand
from numpy import pi, sin, cos, sqrt
    
def jackknife_est(array, func, block_size=1, verbose=False):
    """
    Returns jack_knife est of variance and 1 std err deviation of variance of an array
    Parameters:
        array (numpy array): Array of values to have var estimated
        func: func that takes an array, returns a single value
        block_size (Optional int): Block size of jackknife. default is single value
    Returns:
        var, std_dev_var: variance and 1 std dev of variance as estimated with jackknife.    
    """
    n = np.size(array)
    nb = int(n/block_size)  #number of blocks
    est = func(array)

    jack_pseudos = np.array([nb*est - (nb-1)*func(np.concatenate((array[:i],array[i+block_size:]))) for i in range(nb)])
    jack_est = (1/float(nb))*np.sum(jack_pseudos)  
    jack_est_var = (1/float(nb-1))*np.sum((jack_pseudos-est)**2)
    jack_std_dev_of_est = sqrt(jack_est_var)/sqrt(nb)
    if verbose is True:
        print("Jackknife estimate: {:.3g}".format(jack_est))
        print("1 sig error on Jackknife estimate: {:.3g}".format(jack_std_dev_of_var))
    
    return jack_est, jack_std_dev_of_est
    
def resample_est(array, func, block_size=1, verbose=False):
    """
    Returns jack_knife est of variance and 1 std err deviation of variance of an array
    Parameters:
        array (numpy array): Array of values to have var estimated
        func: func that takes an array, returns a single value
        block_size (Optional int): Block size of jackknife. default is single value
    Returns:
        var, std_dev_var: variance and 1 std dev of variance as estimated with jackknife.    
    """
    n = np.size(array)
    nb = int(n/block_size)  #number of blocks
    est = func(array)

    jack_pseudos = np.array([func(np.concatenate((array[:i],array[i+block_size:]))) for i in range(nb)])
    jack_est = np.mean(jack_pseudos)  
    jack_est_var = np.var(jack_pseudos)
    jack_std_dev_of_est = sqrt(jack_est_var)/sqrt(nb)
    if verbose is True:
        print("Jackknife estimate: {:.3g}".format(jack_est))
        print("1 sig error on Jackknife estimate: {:.3g}".format(jack_std_dev_of_var))
    
    return jack_est, jack_std_dev_of_est

#tests    
if __name__=='__main__':
    samples = 200
    data = rand.normal(loc=0.0, scale=2.0, size=samples)   
    
    jack_std_dev_of_var = jackknife_est(data, np.var, verbose=True)[1]
    frac_diff = jack_std_dev_of_var/( np.var(data)*sqrt(2/(samples-1)) )
    if frac_diff>0.9 and frac_diff<1.1:
        print(frac_diff)
        print("Standard 1 sig error on estimate: {:.3g}".format(np.var(data)*sqrt(2/(samples-1))))
        print("Jacknife 1 std dev of var, test with Gaussian: OK")
    else:
        print(frac_diff)
        print("Standard 1 sig error on estimate: {:.3g}".format(np.var(data)*sqrt(2/(samples-1))))
        print("Jacknife 1 std dev of var, test with Gaussian: FAIL")
    