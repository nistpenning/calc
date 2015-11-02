__author__ = 'jwbritto'

import matplotlib as mpl
import matplotlib.pylab as plt
import sys

# assumes path to source cloned from github.com/nistpenning is in python path
import calc.tools.xkcd_rgb as xrgb


def set_rc(mode='qunat'):
    """ Set matplotlib RC parameters for different plotting applications.

    :param mode: {'default', 'quant', 'seaborn'}
    :return: none
    """

    if mode == 'default':
        plt.rcdefaults()
        return
    elif mode == 'quant':
        style = {'axes.grid': True, 'xtick.major.size': 5,
            'ytick.major.size': 5,
            'savefig.dpi': 80,
            'lines.antialiased': True}
    elif mode == 'clean':
        style = {'axes.grid': False,
            'xtick.major.size': 2,
            'ytick.major.size': 2,
            'savefig.dpi': 80,
            'lines.antialiased': True}
    elif mode == 'seaborn':
        style = {'lines.antialiased': True,
            'text.color': '.15',
            'grid.color': '.8',
            'savefig.dpi': 200,
            'axes.axisbelow': True, 'axes.labelcolor': '.15',
            'xtick.minor.size': 0.0,
            'lines.solid_capstyle': 'round',
            'lines.linewidth': 0.5,
            'axes.edgecolor': '.8', 'font.family': ['sans-serif'],
            'image.cmap': 'Greys', 'ytick.direction': 'out',
            'legend.numpoints': 1,
            'legend.frameon': False,
            'xtick.direction': 'out', 'legend.scatterpoints': 1,
            'axes.linewidth': 1.0, 'ytick.major.size': 0.0,
            'axes.grid': True, 'grid.linestyle': '-',
            'ytick.minor.size': 0.0, 'ytick.color': '.15',
            'font.sans-serif': ['Arial', 'Liberation Sans', 'Bitstream Vera Sans', 'sans-serif'],
            'font.size': 9.0,
            'figure.facecolor': 'white',
            'figure.figsize': (3.37, 2.5),
            'figure.autolayout': True,
            'figure.dpi': 400,
            'xtick.major.size': 0.0,
            'axes.facecolor': 'white', 'xtick.color': '.15'}
    elif mode == 'PRL':
        style = {
            'figure.figsize': (3.375, 2.5), # PRL, column width is 3.375
            'figure.autolayout': True,
            'figure.dpi': 150,
            'figure.facecolor': 'white',
            'figure.subplot.wspace': 0,
            'figure.subplot.hspace': 0,

            'savefig.dpi': 600, # PRL requests 600 dpi

            'font.family': ['serif'],
            'font.size': 10.0,  # PRL, minimum is 7 pt

            'lines.antialiased': True, 'lines.solid_capstyle': 'round',
            'lines.linewidth': 1.0,  # PRL requires >= 0.5 pt

            'axes.grid': True,
            'axes.axisbelow': True,
            'axes.linewidth': 1.0,
            'axes.facecolor': 'white',
            'axes.labelsize': 10,

            'grid.color': '0',
            'grid.linestyle': '-',

            'xtick.minor.size': 2, 'ytick.minor.size': 2,
            'xtick.major.size': 4, 'ytick.major.size': 4,
            'xtick.direction': 'in', 'ytick.direction': 'in',
            'xtick.labelsize': 8, 'ytick.labelsize': 8,

            'image.cmap': 'Greys',

            'legend.numpoints': 1,
            'legend.frameon': False,
            'legend.scatterpoints': 1,
            'legend.fontsize': 8,
            'legend.labelspacing': 0.05,
            'legend.frameon': True,
            'legend.borderaxespad': 0,
            'legend.handlelength': 1.1
        }
    else:
        print('setrc() unknown mode')
        sys.exit(-1)

    for key, value in style.items():
        mpl.rcParams[key] = value

