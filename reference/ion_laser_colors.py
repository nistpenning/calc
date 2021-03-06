__author__ = 'jwbritto'

import sys
import numpy as np
import matplotlib.pylab as plt
import matplotlib
from matplotlib.transforms import Bbox, TransformedBbox, \
     blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector,\
     BboxConnectorPatch

class Laser:
    """Laser properties
    """
    def __init__(self, lambdas, bw, label, type):
        self.type = type
        self.label = label
        self.bw = bw
        if type=="low_power":
            self.l1 = np.array([ [x-bw/2, x+bw/2] for x in lambdas])
        if type=="high_power":
            self.l1 = np.array([lambdas[0]-bw/2, lambdas[0]+bw/2])
            self.l2 = self.l1/2
            self.l3 = self.l1/3
            self.l4 = self.l1/4

    @classmethod
    def high_power(cls, lambda0, bw, label):
        """high power lasers support harmonics

        :param lambda0: center wavelength, float  in nm
        :param bw: gain bandwidth, float  in nm
        :param label: laser label, str
        """
        return cls([lambda0], bw, label, type="high_power")
    @classmethod
    def low_power(cls, lambdas, bw, label):
        """lower power lasers (direct only)

        :param lambdas: center wavelengths, [float]  in nm
        :param bw: tuning bandwidth, float  in nm
        :param label: laser label, str
        """
        return cls(lambdas, bw, label, type="low_power")

    def __str__(self):
        s = "{label}:\n".format(label=self.label)
        if self.type == "low_power:":
            s = s + "  L1={} nm\n".format(self.l0)
            return s
        if self.type == "high_power:":
            s = s + "  L1={} nm\n".format(self.l1)
            s = s + "  L2={} nm\n".format(self.l2)
            s = s + "  L3={} nm\n".format(self.l3)
            s = s + "  L4={} nm\n".format(self.l4)
            return s
        else:
            return ""

class Species:
    """Atomic species wavelengths
    """
    def __init__(self, name, colors):
        """
        :param name: label for ion
        :param colors: dictinary of {{wavelength: "descriptor"}, ...}
        """
        self.name = name
        self.colors = colors

    def __str__(self):
        stmp = ["{}:,{:.1f}nm, ".format(v, k) for k,v in self.colors.items()]
        s = "{name} :: {colors}".format(name=self.name, colors=stmp)
        return s

def connect_bbox(bbox1, bbox2,
                 loc1a, loc2a, loc1b, loc2b,
                 prop_lines, prop_patches=None):
    """

    :param bbox1:
    :param bbox2:
    :param loc1a:
    :param loc2a:
    :param loc1b:
    :param loc2b:
    :param prop_lines:
    :param prop_patches:
    :return:
    """
    #http://matplotlib.org/examples/pylab_examples/axes_zoom_effect.html
    if prop_patches is None:
        prop_patches = prop_lines.copy()
        prop_patches["alpha"] = prop_patches.get("alpha", 1)*0.2

    c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
    c1.set_clip_on(False)
    c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
    c2.set_clip_on(False)

    bbox_patch1 = BboxPatch(bbox1, **prop_patches)
    bbox_patch2 = BboxPatch(bbox2, **prop_patches)

    p = BboxConnectorPatch(bbox1, bbox2,
                           #loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                           loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
                           **prop_patches)
    p.set_clip_on(False)

    return c1, c2, bbox_patch1, bbox_patch2, p


def zoom_effect01(ax1, ax2, xmin, xmax, **kwargs):
    """
    ax1 : the main axes
    ax1 : the zoomed axes
    (xmin,xmax) : the limits of the colored area in both plot axes.

    connect ax1 & ax2. The x-range of (xmin, xmax) in both axes will
    be marked.  The keywords parameters will be used ti create
    patches.
    http://matplotlib.org/examples/pylab_examples/axes_zoom_effect.html
    """

    trans1 = blended_transform_factory(ax1.transData, ax1.transAxes)
    trans2 = blended_transform_factory(ax2.transData, ax2.transAxes)

    bbox = Bbox.from_extents(xmin, 0, xmax, 1)

    mybbox1 = TransformedBbox(bbox, trans1)
    mybbox2 = TransformedBbox(bbox, trans2)

    prop_patches=kwargs.copy()
    prop_patches["ec"]="none"
    prop_patches["alpha"]=0.2

    c1, c2, bbox_patch1, bbox_patch2, p = \
        connect_bbox(mybbox1, mybbox2,
                     loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                     prop_lines=kwargs, prop_patches=prop_patches)

    ax1.add_patch(bbox_patch1)
    ax2.add_patch(bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)

    return c1, c2, bbox_patch1, bbox_patch2, p


def zoom_effect02(ax1, ax2, **kwargs):
    """
    ax1 : the main axes
    ax1 : the zoomed axes

    Similar to zoom_effect01.  The xmin & xmax will be taken from the
    ax1.viewLim.
    http://matplotlib.org/examples/pylab_examples/axes_zoom_effect.html
    """

    tt = ax1.transScale + (ax1.transLimits + ax2.transAxes)
    trans = blended_transform_factory(ax2.transData, tt)

    mybbox1 = ax1.bbox
    mybbox2 = TransformedBbox(ax1.viewLim, trans)

    prop_patches=kwargs.copy()
    prop_patches["ec"]="none"
    prop_patches["alpha"]=0.2

    c1, c2, bbox_patch1, bbox_patch2, p = \
        connect_bbox(mybbox1, mybbox2,
                     loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                     prop_lines=kwargs, prop_patches=prop_patches)

    ax1.add_patch(bbox_patch1)
    ax2.add_patch(bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)

    return c1, c2, bbox_patch1, bbox_patch2, p


def plotit_all(lasers, ions):
    plt.clf()
    fig = plt.figure(1, figsize=(10,2))
    plt.subplots_adjust(left=0.25, right=0.95, top=0.95, bottom=0.15)
    ax1 = plt.subplot(111)

    # plot ion colors
    i=0
    colorlist = ['blue','green','red','cyan','magenta','orange','yellow']
    for ion in ions:
        for nm, label in ion.colors.items():
            s = "{} {}".format(ion.name, label)
            line = ax1.axvline(x=nm, color=colorlist[i], linewidth=2, alpha=0.5)
            #txt = ax1.text(x=nm, y=i+.5, s=s, rotation=90")
        line.set_label(ion.name)
        i+=1
    ax1.legend(loc='upper center',mode='expand', ncol=i)

    # plot laser colors
    i = 0
    xmax = 0
    yticks = []; yticklabels = []
    xticks = []; xticklabels = []
    laser_alpha = 1
    for r in lasers:
        yticks.append(0.5 + i)
        yticklabels.append(r.label)
        # fundamental
        if r.type == "low_power":
            # there's a list of wavelengths for a low power laser class
            xbars = [(x[0], x[1]-x[0]) for x in r.l1]
            xmax = np.max([xmax, np.max(r.l1)])
            ax1.broken_barh(xbars , (i, 1), facecolors='grey',
                        alpha=laser_alpha, linewidth=0)
        if r.type == "high_power":
            xbars = [(r.l1[0], r.l1[1]-r.l1[0])]
            xmax = np.max([xmax, r.l1[1]])
            b1 = ax1.broken_barh(xbars , (i, 1), facecolors='red',
                        alpha=laser_alpha, linewidth=0,
                            label="fundamental")

        # second and higher harmonics only for high-power lasers
        if r.type == "high_power":
            # second harmonic
            xbars = [(r.l2[0], r.l2[1]-r.l2[0])]
            b2 = ax1.broken_barh(xbars, (i, 1), facecolors='green',
                            alpha=laser_alpha, linewidth=0,
                            label="2nd harmonic")
            # third harmonic
            xbars = [(r.l3[0], r.l3[1]-r.l3[0])]
            b3= ax1.broken_barh(xbars, (i, 1), facecolors='blue',
                            alpha=laser_alpha, linewidth=0,
                            label="3rd harmonic")
            # fourth harmonic
            xbars = [(r.l4[0], r.l4[1]-r.l4[0])]
            b4 = ax1.broken_barh(xbars, (i, 1), facecolors='purple',
                            alpha=laser_alpha, linewidth=0,
                            label='4th harmonic')
        i += 1

    yticks.append(1 + i)
    yticklabels.append("Sources:")

    # spruce up axes
    ax1.grid(True)
    ax1.set_yticks(yticks)
    ax1.set_yticklabels(yticklabels)
    ax1.set_xlim([125, xmax])
    ax1.set_xlabel("wavelength [nm]")
    plt.show()


def plotit_uv_zoom(lasers, species, fig_title=''):
    nspecies = len(species)
    nions = len(ion_qubits)

    xstretch = 0.1
    zoom_x_max = 425
    xmin = 225
    plt.clf()
    plt.figure(1, figsize=(5,5))
    plt.subplots_adjust(left=0.25, right=0.95, top=0.70, bottom=0.1)
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)

    i = 0
    yticks = []; yticklabels = []
    xticks = []; xticklabels = []

    xmax = 0
    for r in lasers:
        yticks.append(0.5 + i)
        yticklabels.append(r.label)
        # fundamental
        if r.type == "low_power":
            # there's a list of wavelengths for a low power laser class
            xbars = [(x[0], x[1]-x[0]) for x in r.l1]
            xmax = np.max([xmax, np.max(r.l1)])
        if r.type == "high_power":
            xbars = [(r.l1[0], r.l1[1]-r.l1[0])]
            xmax = np.max([xmax, r.l1[1]])
        ax1.broken_barh(xbars , (i, 1), facecolors='red')
        ax2.broken_barh(xbars , (i, 1), facecolors='red')

        # second and higher harmonics only for high-power lasers
        if r.type == "high_power":
            # second harmonic
            xbars = [(r.l2[0], r.l2[1]-r.l2[0])]
            ax1.broken_barh(xbars , (i, 1), facecolors='green')
            ax2.broken_barh(xbars , (i, 1), facecolors='green')
            # third harmonic
            xbars = [(r.l3[0], r.l3[1]-r.l3[0])]
            ax1.broken_barh(xbars , (i, 1), facecolors='blue')
            ax2.broken_barh(xbars , (i, 1), facecolors='blue')
            # fourth harmonic
            xbars = [(r.l4[0], r.l4[1]-r.l4[0])]
            ax1.broken_barh(xbars , (i, 1), facecolors='purple')
            ax2.broken_barh(xbars , (i, 1), facecolors='purple')
        i += 1

    # plot ion colors
    lbl_sqz = 3 # labels are too close, nm
    sqzd = False
    all_colors = []  # list of all ion colors
    for ion in species:
        for nm, label in ion.colors.items():
            sqzd = any( [(abs(x-nm)<lbl_sqz) or (abs(nm-x)<lbl_sqz)
                         for x in all_colors] )
            all_colors.append(nm)
            if sqzd == False:
                ys = nions-0.75
                ynm = -0.75
            else:
                ys = nions+1.7
                ynm = .25
                sqzd = False
            s = "{}: {}".format(ion.name, label)
            ax1.axvline(x=nm, color='y')
            ax1.text(x=nm, y=ys, s=s, rotation=90, size=9,
                     verticalalignment='bottom', horizontalalignment='center')
            ax1.text(x=nm, y=ynm, s=nm, rotation=90, size=9,
                     verticalalignment='bottom', horizontalalignment='center')

    # get axis labels sorted out
    if nspecies:
        ax1.xaxis.set_tick_params(labeltop='on')
        ax1.xaxis.set_tick_params(labelbottom='off')
        ax1.set_xticks(xticks)
        ax1.set_xticklabels(xticklabels)
        for label in ax1.xaxis.get_ticklabels():
            label.set_rotation(90)
        caption_s = "LC is laser cooling\n" \
            "RP is repump\n" \
            "RA is Raman-transition\n" \
            "PI is photo-ionization"
        plt.figtext(x=0.01, y=0.99, s=caption_s,
            verticalalignment='top')

    else:
       caption_right="RED is fundamental\n" \
                  "   (>500mW tuning)\n" \
                  "GREEN is doubled\n" \
                  "BLUE is trippled\n" \
                  "PURPLE is quadrupled"
       plt.figtext(x=0.75, y=0.99, s=caption_right,
            horizontalalignment="left", verticalalignment='top')
    ax1.grid(True)
    ax1.set_yticks(yticks)
    ax1.set_yticklabels(yticklabels)
    ax2.set_xlim([xmin*(1-xstretch), xmax*(1+xstretch)])
    ax2.set_ylim([0,i+0.5])
    ax2.grid(True)
    ax2.set_yticks(yticks)
    ax2.set_yticklabels(yticklabels)
    ax2.set_ylim([0,i+0.5])
    ax2.set_xlabel("wavelength [nm]")
    ax1.set_xlim([xmin, zoom_x_max])
    zoom_effect01(ax1, ax2, xmin, zoom_x_max)
    plt.figtext(x=.5, y=0.99, s=fig_title,
                horizontalalignment="center", verticalalignment='top',
                weight='bold')
    plt.show()

# diode source data from
# http://tf.boulder.nist.gov/general/pdf/2765.pdf
#
lasers_nist = [Laser.high_power(1118, 40, "1118nm OPSL"),
        Laser.high_power(1156, 60, "1156nm OPSL K3381"),
        Laser.high_power(1200, 40, "1200nm OPSL"),
        Laser.high_power(705, 30, "705nm OPSL wish"),
        Laser.low_power([0],
                        2, ""),
        # Laser.low_power([1083, 671, 640, 766, 850, 854, 866, 649, 658, 812,
        #                  780, 795, 1033, 1092, 882, 650, 935, 718],
        #                 2, "direct sources"),
        # Laser.low_power(np.array([626, 1178, 844, 852, 922, 814, 842, 1108, 910, 986])/2,
        #                 2, "doubled sources"),
        # Laser.low_power(np.array([939, 855, 840, 1179, 1191, 984, 1197, 984, 1107])/3,
        #                 2, "tripled sources"),
        # Laser.low_power(np.array([940, 1140, 1120, 916, 856, 908])/4,
        #                 2, "quadrupled sources")
        ]
[print(x) for x in lasers_nist]

# import OPSL chips
# not for public distribution
opsl_path = "C:\\Users\\jwbritto\\Google Drive\\0workSync\\reference\\tomi_laser_colors_opsl.py"
exec(open(opsl_path).read())

# qubit and clock ions
# from ref [0] unless otherwise noted
# [0] http://tf.boulder.nist.gov/general/pdf/2765.pdf
# [1] http://link.aps.org/doi/10.1103/PhysRevA.85.012502
ion_qubits = [Species("Be+", {313:"LC,RA",
                               235:"PI"}),
        Species("Mg+", {280:"LC,RA",
                        285:"PI"}),
        Species("Yb+", {297:"PI,LC[1]",
                        399:"PI",
                        556:"PI",
                        328:"LC,RA",
                        369:"LC,RA",
                        935:"RP"}),
        Species("Ca+", {422:"PI",
                        850:"RP",
                        854:"RP",
                        866:"RP",
                        393:"LC,RA",
                        397:"LC,RA"}),
        Species("Sr+", {461:"PI",
                        407:"LC,RA",
                        421:"LC,RA",
                        1033:"RP",
                        1092:"RP"}),
        Species("Ba+", {554:"PI",
                        455:"LC,RA",
                        493:"LC,RA",
                        650:"RP"})]
[print(x) for x in ion_qubits]

species_empty = []

neutrals = [Species("He*", {1083:"LC,DP"}),
               Species("Li", {671:"LC,RA,DP"})]

#plotit_all(lasers, ions)
# plotit_uv_zoom(lasers_nist, ion_qubits,
#                fig_title="OPSLs at NIST (Burd)")
# plotit_uv_zoom(lasers_tomi_tested, species_empty,
#                fig_title="OPSLs Demonstrated by tut.fi (2015)")
plotit_uv_zoom(lasers_tomi_prospects, species_empty,
               fig_title="Anticipated OPSLs (2015, tut.fi)")

