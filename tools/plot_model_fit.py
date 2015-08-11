# -*- coding: utf-8 -*-
"""
Created on Tue Jun 02 10:19:20 2015

@author: jgb
"""

import numpy as np
from numpy import pi, sin, cos
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.text
import scipy.optimize as opt

def auto_extent(x,y):
    xm = np.min(x)
    ym = np.min(y)
    xx = np.max(x)
    yx = np.max(y)
    full_x = abs(xx - xm)
    full_y = abs(yx - ym)
    axis = [xm-0.05*full_x, xx + 0.1*full_x, ym-0.1*full_y, yx+0.1*full_y]
    return axis

def save_data_txt(filename, out_list):
    """
    filename: string with filename, including extenstion
    out_array: list with different data to be output
    """
    s=''
    fh = open(filename, 'w+')
    out = map(list, zip(*out_list))
    for row in out:
        for i in row:
            s = s+'%f,'% (i)
        s = s+" \n"
    fh.write(s)
    fh.close()

def plot_fit(x,y,fitfunc,fitguess,
            yerr=None,
            hold=[],
            labels=['X','Y','default'],
            axis=None,
            save=False,
            show=True,
            fmt_data='o',
            fmt_fit='-'):
    """Plot & fit with supplied model.

    :param x: numpy.array
    :param y: numpy.array
    :param fitfunc: fit function for N free parameters f(x, k0, k1, ..., kN); returns a scalar
    :param fitguess: guess for each kN [guess_k0, guess_k1, ..., guess_kN]
    :param yerr: error for y [numpy.array]
    :param hold: model parameters to hold fixed [b1, b2, ..., bN]; 'none'; 'all'
    :param labels: ["xlabel", "ylabel", "plot title"]
    :param axis: axis extent [xmin, xmax, ymin, ymax]; 'auto'
    :param save: save as .png, True or False
    :param show: show data and fit on plot, True or False
    :param fmt_data: matplotlib plot format string; default 'o'
    :param fmt_fit: matplotlib plot format string; default '-'
    :return: [[fit_param], [fit_error]]
        fit_param is vector of fit coefficients of length N
        fit_error sqrt of the diagonals of the covariance matrix (1 sigma confidence interval)
    """

    #process hold parameter to make a fit model with specified number of free params
    if hold=='all':
        hold = np.ones(np.size(fitguess), dtype=bool)

    if np.size(hold) == 0:
        #default, fit it all
        hold = np.zeros(np.size(fitguess), dtype=bool)
        varin = fitguess
    else:
        #create the list of free parameter guesses
        #note ~ performs the NOT function on the boolean array
        varin = fitguess[~hold]

    if hold.all():
        #No free params, just plot the model, don't perform fit
        print("No free parameters, plotting model")
        x_curve = np.linspace(np.min(x),np.max(x), num=200)
        curve_fit = fitfunc(x_curve, *fitguess)
        popt = fitguess
        perr = np.zeros(np.size(popt))
    else:
        #dynamically define the fit function with subset of fixed arguments
        def func(x,*var):
            args = np.array([])
            var_count = 0
            #build array of parameters, based on the hold paramter
            for i in range(np.size(fitguess)):
                if hold[i]:
                    #get from fitguess, paramter not varied
                    args = np.append(args,fitguess[i])
                else:
                    #get from input to func, will actaully be varied
                    args = np.append(args,var[var_count])
                    var_count+=1
            return fitfunc(x,*args)
        #actually do the fit, with the subset of params, varin
        popt,pcov = opt.curve_fit(func, x, y, p0=varin, sigma=yerr)
        try:
            perr = np.sqrt(np.diag(pcov))
        except ValueError:
            print("Ill defined fits")
            perr = np.zeros(np.size(fitguess))
        x_curve = np.linspace(np.min(x),np.max(x), num=200)
        curve_fit = func(x_curve, *popt)

    poptf = ('%.3g, ' * len(popt))[:-1] % tuple(popt)
    perrf = ('%.3g, ' * len(perr))[:-1] % tuple(perr)
    fit_message = "fit: {} +/- {} (1 sigma)".format(poptf, perrf)

    if show is True:
        # build figure
        if axis == None:
            axis = [0.0, 1.1*np.max(x), 0.0,  1.1*np.max(y)]
        elif axis == 'auto':
            axis = auto_extent(x,y)
        plt.axis(axis)
        plt.errorbar(x,y,yerr=yerr,fmt=fmt_data)
        plt.plot(x_curve,curve_fit,fmt_fit)

        # labels
        plt.xlabel(labels[0])
        plt.ylabel(labels[1])

        # fit message on plot
        ym = axis[-2]
        xm = axis[0]
        y_pos_msg = ym-(0.4*np.abs(axis[3]-ym))
        if xm is 0.0:
            x_pos_msg = 0.0
        else:
            x_pos_msg = xm
        plt.figtext(x_pos_msg, y_pos_msg, fit_message, fontsize=9)

        name = fitfunc.__name__
        if save:
            name_out = name+'_fit.png'
            plt.savefig(name_out, format='png', dpi=300,
                        transparent=True, bbox='tight')

        if labels[2]=='default':
            plt.title(name)
        else:
            plt.title(labels[2])

    else:
        print(fit_message)

    return popt, perr

def plot_polyfit(x,y, fitguess,
            yerr=None,
            hold=[],
            labels=['X','Y','default'],
            axis='default',
            save=False,
            show=True,):
    """Plot data and the fit to a polynominal, with ability to fix parameters.

    Parameters
    ----------
    x : array
        x data for plot
    y : array
        y data for plot        
    fitguess : array
        polynominal order to fit to k0+x*k1+x^2*k2...
        the number of guesses supplied will determine polynominal order

    Keyword arguments
    -----------------
    yerr : array
        1 std deviation for the y data to weight the fits
        default is equal weight
    hold : list
        list of booleans, what model parameters to hold fixed (default none)
        can also pass hold='all' to hold all params fixed
    labels : list
        list of strings, xlabel ylabel, plot title
    axis : list
        list of axis extent, [xmin, xmax, ymin, ymax] (default 0, xmax, 0, ymax)
    save: boolean
        boolean to choose if figure is saved as png file (dafault no save)
    show: boolean
        boolean to choose if the data and fit are shown on a plot
    Returns
    -------
    popt : array
        the fitted values only (fitguess if no free parameters)
    perr : array
        sqrt of the diagonals of the covariance matrix (1 sigma confidence interval)

    """
# replace the old fitfunc with polynominal fit
    def fitfunc(x, *params):
        return sum([p*(x**i) for i, p in enumerate(params)])

#process hold parameter to make a fit model with specified number of free params
    if hold=='all':
        hold = np.ones(np.size(fitguess), dtype=bool)

    if np.size(hold) == 0:
        #default, fit it all
        hold = np.zeros(np.size(fitguess), dtype=bool)
        varin = fitguess
    else:
        #create the list of free parameter guesses
        #note ~ performs the NOT function on the boolean array
        varin = fitguess[~hold]

    if hold.all():
        #No free params, just plot the model, don't perform fit
        print("No free parameters, plotting model")
        x_curve = np.linspace(np.min(x),np.max(x), num=200)
        curve_fit = fitfunc(x_curve, *fitguess)
        popt = fitguess
        perr = np.zeros(np.size(popt))
    else:
        #dynamically define the fit function with subset of fixed arguments
        def func(x,*var):
            args = np.array([])
            var_count = 0
            #build array of parameters, based on the hold paramter
            for i in range(np.size(fitguess)):
                if hold[i]:
                    #get from fitguess, paramter not varied
                    args = np.append(args,fitguess[i])
                else:
                    #get from input to func, will actaully be varied
                    args = np.append(args,var[var_count])
                    var_count+=1
            return fitfunc(x,*args)
        #actually do the fit, with the subset of params, varin
        popt,pcov = opt.curve_fit(func, x, y, p0=varin, sigma=yerr)
        try:
            perr = np.sqrt(np.diag(pcov))
        except ValueError:
            print("Ill defined fits")
            perr = np.zeros(np.size(fitguess))
        x_curve = np.linspace(np.min(x),np.max(x), num=200)
        curve_fit = func(x_curve, *popt)

#fit message
    poptf = ('%.3g, ' * len(popt))[:-1] % tuple(popt)
    perrf = ('%.3g, ' * len(perr))[:-1] % tuple(perr)

    fit_message = 'Curve fit results: ' + poptf
    fit_mess2 = '\n    uncertianties: ' +  perrf

    fit_message = fit_message + fit_mess2

    if show is True:
        #build figure
        plt.close()

        if axis == 'default':
            axis = [0.0, 1.1*np.max(x), 0.0,  1.1*np.max(y)]
        elif axis == 'auto':
            axis = auto_extent(x,y)
        plt.axis(axis)

        if yerr is None:
            plt.plot(x,y,'o')
        else:
            plt.errorbar(x,y,yerr=yerr,fmt='o')
        plt.plot(x_curve,curve_fit,'-')

        #labels
        plt.xlabel(labels[0])
        plt.ylabel(labels[1])
        ym = axis[-2]
        xm = axis[0]
        y_pos_msg = ym-(0.3*np.abs(axis[3]-ym))

        if xm is 0.0: x_pos_msg = 0.0
        else: x_pos_msg = xm

        plt.text(x_pos_msg, y_pos_msg, fit_message, fontsize=10)

        name = fitfunc.__name__
        if save:
            name_out = name+'_fit.png'
            plt.savefig(name_out, format='png', bbox='tight')

        if labels[2]=='default':
            plt.title(name)
        else:
            plt.title(labels[2])

        plt.show()
    else:
        print(fit_message)

    return popt, perr
