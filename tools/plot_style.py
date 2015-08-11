# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 20:22:24 2015

@author: jgb
"""

# define some plot colors to match my matplotlibrc
red = '#A60628'
blue = '#348ABD'
purple = '#7A68A6'
green = '#467821'
orange = '#D55E00'
pink = '#CC79A7'
cyan = '#56B4E9'
aqua = '#009E73'
yellow = '#F0E442'
navy = '#002b36'

def save_data_txt(filename, out_list, col_names=False):
    """
    filename: string with filename, including extenstion
    out_array: list with different data to be output
    col_names: string with comma seperated column names
    """
    cn = col_names + " \n"
    fh = open(filename, 'w+')    
    if col_names is not False:
        fh.write(cn)

    s=''
    out = map(list, zip(*out_list))
    for row in out:
        for i in row:
            s = s+'%f,'% (i)
        s = s+" \n"
    fh.write(s)
    fh.close()