# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 21:08:35 2015

@author: justinbohnet
"""

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

def plot_fit(x,y,fitfunc,fitguess,
            hold=[],
            labels=['X','Y','default'],
            axis='default',
            save=False):
    """Plot data and the fit to a given model, with ability to fix parameters.

    Parameters
    ----------
    x : array
        x data for plot
    y : array
        y data for plot
    fitfunc : function
        function of the form f(x,*args)
    fitguess : array
        must supply values for all the parameters in the model given by fitfunc

    Keyword arguments
    -----------------
    hold : list
        list of booleans, what model parameters to hold fixed (default none)
    labels : list
        list of strings, xlabel ylabel, plot title
    axis : list
        list of axis extent, [xmin, xmax, ymin, ymax] (default 0, xmax, 0, ymax)
    save: boolean
        boolean to choose if figure is saved as png file (dafault no save)
    Returns
    -------
    popt : array
        the fitted values only (fitguess if no free parameters)
    perr : array
        sqrt of the diagonals of the covariance matrix (1 sigma confidence interval)

    """


#process hold parameter to make a fit model with specified number of free params
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
        popt,pcov = opt.curve_fit(func, x, y, p0=varin)
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

#build figure
    plt.close()

    if axis == 'default':
        axis = [0.0, 1.1*np.max(x), 0.0,  1.1*np.max(y)]
    elif axis == 'auto':
        xm = np.min(x)
        ym = np.min(y)
        full_x = abs(np.max(x) - xm)
        full_y = abs(np.max(y) -ym)
        axis = [xm-0.1*full_x, 1.1*np.max(x), ym-0.1*full_y, 1.1*np.max(y)]
    plt.axis(axis)

    plt.plot(x,y,'o')
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

    return popt, perr