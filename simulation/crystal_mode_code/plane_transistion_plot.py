__author__ = 'sbt'
"""
Makes a plot of the rotation frequency of the
2-1 plane transistion for a given configuration of the Ion trap.
"""
from mode_analysis_code import ModeAnalysis
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

if __name__ == "__main__":

    # Select the trapping and wall potentials which will be used
    # for all future runs
    trappotential = (0.0, -873, -1000)
    wallpotential = 0.1
    precision_solving = True

    # Determines the number of ions to find the transition frequency for.
    nionlist = [19, 20, 26, 37, 50, 61, 75, 91, 110, 127, 130, 169, 190, 217, 231, 300, 331]

    currentfrequency = 93

    transistionfrequencies = []

    # Iterate through number of ions to test for stability
    for N in nionlist:

        if N > 100:
            # Set to false to decrease run time for the biggest crystals
            # and run with un-perturbed crystals
            # (potentially not global energy minimum)
            precision_solving = True

        # Instantiate a crystal and see if it is stable
        crystal = ModeAnalysis(N=N, Vtrap=trappotential, Ctrap=1.0, ionmass=None,
                               B=4.4588, frot=currentfrequency, Vwall=wallpotential,
                               wall_order=2, quiet=False, precision_solving=precision_solving)
        crystal.run()

        # Increase the frequency until stability is lost- most important for the first
        # crystal tested
        while crystal.is_plane_stable():
            print("Crystal of", N, "is currently at", currentfrequency,
                  "increasing to ", currentfrequency + 1)

            currentfrequency += 1
            crystal = ModeAnalysis(N=N, Vtrap=trappotential, Ctrap=1.0, ionmass=None,
                                   B=4.4588, frot=currentfrequency, Vwall=wallpotential,
                                   wall_order=2,
                                   quiet=False, precision_solving=precision_solving)

            crystal.run()

        # When frequency is lost, reduce to find when it resumes
        while not crystal.is_plane_stable():
            print("Found turning point: reducing frequency from", currentfrequency, "to ",
                  currentfrequency - 1)
            currentfrequency -= 1

            crystal = ModeAnalysis(N=N, Vtrap=trappotential, Ctrap=1.0, ionmass=None,
                                   B=4.4588, frot=currentfrequency, Vwall=wallpotential,
                                   wall_order=2,
                                   quiet=False, precision_solving=precision_solving)
            crystal.run()

        # Once stability has resumed the lowest frequency at which 1->2 transition occurs is stored
        print("Transistion frequency is", currentfrequency + 1, " for number of ions", crystal.Nion)
        transistionfrequencies.append(currentfrequency + 1)

    print("Transitions found:")
    print("nions:", nionlist)
    print("frequencies", transistionfrequencies)

    #########################################

    transfreq=transistionfrequencies
    nions=nionlist

    shells=[1,2,3,4,5,6,7,8,9,10]
    shelln=[7,19,37,61,91,127,169,217,271,331]

    def func(x, a, b, c):
        return a * np.exp(-b * x) + c

    fig = plt.figure(figsize=(14, 12))
    plt.rcParams['font.size'] = 16
    ax = fig.add_subplot(1,1,1)
    for i in range(len(transfreq)):
        if nions[i] in shelln:
             plt.plot(transfreq[i],nions[i],"o",color='red')
        else:
             plt.plot(transfreq[i],nions[i],"o",color='blue')


    plt.title("1-2 Plane Transistion for $V_{Mid}=-.873, \ V_{Center}=-1.0 \ (kV) V_{Wall} =1 V$", y=1.02)
    plt.xlabel("Transistion Frequency (kHz)")
    plt.ylabel("Number of Ions")


    major_ticks = np.arange(min(transfreq),max(transfreq),2)
    minor_ticks = np.arange(min(transfreq),max(transfreq),.5)
    print(major_ticks)
    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks, minor=True)
    yticks=np.arange(0,400,25)
    yticksmin=np.arange(0,400,5)
    ax.set_yticks(yticks)
    ax.set_yticks(yticksmin, minor=True)

    fig = plt.grid(True)

    fig = plt.xlim([min(transfreq)*.99,max(transfreq)*1.01])

    popt, pcov = curve_fit(func, transfreq, nions,p0=[127,.1,122])
    print(popt)
    x=np.linspace(min(transfreq)*.99,max(transfreq)*1.01,200)
    plt.plot(x, func(x, *popt), 'r-', label="Fitted Curve",color="black")
    plt.legend(loc=1)


    for N in shelln:
        plt.plot([min(transfreq)*.99,max(transfreq)*1.01],[N,N],"--",color='black')
    for N in shells:
        plt.text(max(transfreq)*1.013,shelln[N-1],"%d" %N)
    plt.show()

