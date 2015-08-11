# -*- coding: utf-8 -*-
"""
Created on Wed May 27 10:40:23 2015

@author: justinbohnet
"""

import os
import numpy as np
from numpy import sin, cos, pi, sqrt, exp
import matplotlib.pyplot as plt
#import scipy.stats.norm as norm

import ODF
import hfGUIdata
import plot_model_fit as pt

def OAT(psi, chi, N, t):
    """
    psi: quadrature angle [radians] -- could be array
    params:
    chi:  2* (Jbar) / N
    N: atom number
    t: total interaction time
    """
    Acoef = 1-cos(2*chi*t)**(N-2)
    Bcoef = 4*sin(chi*t)*cos(chi*t)**(N-2)
    delt = 0.5*np.arctan2(Bcoef,Acoef)
    varJz = N/4.*(1+(N/4.-0.25)*(Acoef - sqrt(Acoef**2+Bcoef**2)*cos(2*(psi+delt))))
    return varJz

def sq_analysis(max_c, min_c, N, N_err, sigA, k0, Jbar_1kHz):
    # Get data
    file_name, scandata, m, pmterr, trials, data = hfGUIdata.get_raw_counts()

    #data = np.genfromtxt(file_name,delimiter=",",names=True,dtype=None)
    arm_time = np.mean(data['arm_t'])
    pi_time = np.mean(data['middle_t'])
    det_n = np.mean(data['det_n'])

    # Calculate derived quantities
    detune = det_n*1e3/(arm_time)  # kHz
    phi = (scandata/pi_time*180.0) #degrees
    chi = 2*Jbar_1kHz/detune/N  # interaction strength, scaled for detuning

    pmterr_uncert = pmterr/sqrt(2*trials)
    m_avg = np.mean(m)
    K = max_c - min_c

    def noise_model(phi, Ntot, sigA, k0):
        phi = phi*pi/180
        varPN = K**2 * (sin(phi)**2)/4.0/Ntot  # extra factor of K is for conversion to photons
        varAN = (sigA * sin(phi))**2  # sigA is the std dev of rotation of theta in radians
        varSN = K*(sin(phi/2.)**2) + min_c
        varTN = k0**2 * (K*(sin(phi/2.)**2) + min_c)**2
        return sqrt(varPN + varAN + varSN + varTN)

    # to predict the SQL, don't include the added noise in the anti-squeezed quad
    SQL_stdev = noise_model(90., N, 0.0, k0)
    SQL_plus = noise_model(90., N+N_err, 0.0, k0) * np.ones(np.size(phi))
    SQL_minus = noise_model(90., N-N_err, 0.0, k0) * np.ones(np.size(phi))

    title_name = "data set: "+file_name[13:-8]
    title = "%s, $t_{a}$:%d us, $\delta$:%.3g kHz"%(title_name, arm_time, detune)

    label = ['Final rotation [degrees]', 'Std. Dev. PMT counts', title]

    plt.close()
    plt.errorbar(phi, pmterr, yerr=pmterr_uncert, fmt='o')
    plt.fill_between(phi, SQL_plus, y2=SQL_minus, facecolor=pt.red, alpha=0.5 )
    plt.xlabel(label[0])
    plt.ylabel(label[1])
    plt.title(title)
    plt.axis([0,np.max(phi),0,1.5*np.max(pmterr)])

    # try to calculated the squeezed state from OAT
    psi = np.linspace(0,-720,num=400)*pi/180.
    t = 2*arm_time*1e-6
    """
    Acoef = 1-cos(2*J*t)**(N-2)
    Bcoef = 4*sin(J*t)*cos(J*t)**(N-2)
    delt = 0.5*np.arctan2(Bcoef,Acoef)
    varJz = N/4.0*(1+(N/4.0-0.25)*(Acoef - sqrt(Acoef**2+Bcoef**2)*cos(2*(psi+delt))))
    """
    varJz = OAT(psi, chi, N, t)

    sig_squeeze = sqrt(varJz/N**2 * K**2)
    sig_squeeze = sqrt(sig_squeeze**2 + 0.0**2 + np.mean(m) + (k0*np.mean(m))**2)
    plt.plot(180/pi*(-psi),sig_squeeze,'-',c=pt.blue)
    #plt.axis([-5,370,0,80])
    plt.show()

    # Fitting to estimate of anti-squeeze, angles
    plt.close()
    def varfit(psi, maxf, minf, angle):
        psi = pi*psi/180.0
        return (maxf - minf)*cos(psi + angle)**2 + minf
    res, res_err = pt.plot_fit(phi, pmterr**2, varfit, np.array([4.0e4, 2.5e3, 0.0]),
                               hold=np.array([False,False,False]), axis='auto', show=True)

    # check that a cosine is a valid model for extracting parameters from data
#    plt.plot(180/pi*(psi),sig_squeeze**2,'-')
#    plt.plot(phi, pmterr**2,'o')
#    plt.axis([0,370,0,60e3])

    alpha_deg = res[2]*180.0/pi
    alpha_deg_err = res_err[2]*180.0/pi

    if res[0] > res[1]:
        var_max_est = res[0]
        var_max_err = res_err[0]

        var_min_est = res[1]
        var_min_err = res_err[1]
    else:
        var_max_est = res[1]
        var_max_err = res_err[1]

        var_min_est = res[0]
        var_min_err = res_err[0]

    sdev_min = np.min(pmterr)
    m_min = m[pmterr == sdev_min][0]

    #R = (sdev_min**2 - ((sigA)**2 + m_min + (k0*m_min)**2))/(K**2/4.0/N)
    R = (sdev_min**2 - m_min)/(K**2/4.0/N)  # conservative inferred squeezing
    RO = (sdev_min**2)/(K**2/4.0/N)
    RO_all = 10*np.log10(pmterr**2/(K**2/4.0/N))
    print(pmterr)

    print("=========== Numbers from the data ==============")
    print(title_name)
    print('Average PMT counts: {:.4g}, SN: {:.4g}'.format(m_avg,sqrt(m_avg)))
    print('N: {}, PN: {:.4g} counts'.format(N, K/2./sqrt(N)))
    print('Number of photons for N ions, K: {:.4g}'.format(K))
    print('Predicted SQL (PN+SN+TN): {:.4g}'.format(SQL_stdev))
    print('predicted chi from ACSS: {}'.format(J))
    #print('Error bars are statistical: 1/sqrt(2*trials)')
    print('---------------------------------------------------')
    print("Lowest Std Dev. {:.4g} counts".format(sdev_min))
    print("Observed spin noise reduction: {:.4g}, {:.4g}".format(RO, 10*np.log10(RO)))
    print("Inferred spin noise reduction (subtract SN only): {:.4g}, {:.4g}".format(R, 10*np.log10(R)))
    print('Angle alpha [deg]: {0:.4g}'.format(alpha_deg) +  '+-' + '{0:.4g}'.format(alpha_deg_err))
    print('Anti-squeezed std dev. {0:.4g} +- {1:.4g} counts'.format(sqrt(var_max_est),sqrt(var_max_err)))

    return 2*arm_time, detune, m, np.min(pmterr), phi, pmterr, psi, sig_squeeze, RO

def cal_analysis(mask_range, hold=False, Nguess=100, show_var=True, show_avg=False):

    # Get data
    file_name, scandata, m, pmterr, trials, data = hfGUIdata.get_raw_counts()

    #data = np.genfromtxt(file_name,delimiter=",",names=True,dtype=None)

    pi_time = hfGUIdata.get_ionProp_value('sf%fitParams%sf_fitParam_tpi')
    det_t = hfGUIdata.get_ionProp_value('detection%det_t')

    # Calculate derived quantities
    phi = ((scandata/pi_time) * 180.0)-180.0  # degrees

    # Make plot for calibrating Bloch Vector angle
    title_name = "data set: " + file_name[13:-8]
    title = "%s"%(title_name)
    label = ['Polar angle $\phi$ [degrees]', 'Avg PMT counts', title]

    def constrast_fit(phi, max_counts, min_counts):
        phi = phi*pi/180
        return (max_counts - min_counts)*sin(phi/2.)**2 + min_counts
    res, res_err = pt.plot_fit(phi, m, constrast_fit, np.array([1000.0, 150.0]),
                               labels=label, axis='default', show=show_avg)
    # Store results for PMT counts calibration
    max_c_fit = res[0]
    min_c_fit = res[1]
    K = max_c_fit - min_c_fit  # photons per N ions
    N = Nguess

    #option to make a mask for noise calibration
    if mask_range is 0:
        phi_m = phi
        pmterr_m = pmterr
    else:
        mask = np.ones(np.size(phi)).astype(bool)
        mask[mask_range[0]:mask_range[1]] = False
        phi_m = phi[mask]
        pmterr_m = pmterr[mask]

    label = ['Polar angle $\phi$ [degrees]', 'Std. Dev. PMT counts', title]
    fitguess = np.array([N, 0.01, 0.01])

    def noise_model(phi, Ntot, sigA, k0):
        phi = phi*pi/180
        varPN = K**2 * (sin(phi)**2)/4.0/Ntot  # extra factor of K is for conversion to photons
        varAN = (sigA * sin(phi))**2  # sigA is the std dev of rotation of theta in radians
        varSN = K*(sin(phi/2.)**2) + min_c_fit
        varTN = k0**2 * (K*(sin(phi/2.)**2) + min_c_fit)**2
        return sqrt(varPN + varAN + varSN + varTN)


    if hold is False:
        res, res_err = pt.plot_fit(phi_m, pmterr_m, noise_model, fitguess,
                               labels=label, axis='default', show=False)
    else:
        res, res_err = pt.plot_fit(phi_m, pmterr_m, noise_model, fitguess,
                               labels=label, hold=hold, axis='default', show=False)
    if np.size(res) is 3:
        Nfit = res[0]
        sigA = res[1]
        k0 = res[2]
    if np.size(res) is 2:
        Nfit = N
        sigA = res[0]
        k0 = res[1]

    phi_show = np.linspace(0,360,num=200)
    phi_rad = phi_show*pi/180.0
    varPN = K**2 * (sin(phi_rad)**2)/4.0/Nfit
    varAN = (sigA * sin(phi_rad))**2
    varSN = K*(sin(phi_rad/2.)**2) + min_c_fit
    varTN = k0**2 * (K*(sin(phi_rad/2.)**2) + min_c_fit)**2

    if show_var is True:
        # Show the fitted noise on the plot
        plt.plot(phi_m, pmterr_m**2,'o')
        plt.plot(phi_show, noise_model(phi_show, Nfit, sigA, k0)**2)
        plt.plot(phi_show, (varPN), label='PN')
        plt.plot(phi_show, (varAN), label='AN')
        plt.plot(phi_show, (varSN), label='SN')
        plt.plot(phi_show, (varTN), label='TN')
        plt.xlabel('Polar angle  $\phi$ [degrees]')
        plt.ylabel('Var PMT counts')
        plt.title(title)
        plt.legend()
        plt.show()
    else: pass

    return max_c_fit, min_c_fit, Nfit, sigA, k0

def con_analysis(max_full=100.0, min_full=0.0, show=True):
    print('##### Do contrast analysis ######')
    # Get data
    file_name, scandata, m, pmterr, trials, data = hfGUIdata.get_raw_counts()

    #data = np.genfromtxt(file_name,delimiter=",",names=True,dtype=None)

    pi_time = hfGUIdata.get_ionProp_value('sf%fitParams%sf_fitParam_tpi')
    det_t = hfGUIdata.get_ionProp_value('detection%det_t')

    # Calculate derived quantities
    phi = ((scandata/pi) * 180.0)  # degrees

    # Make plot for calibrating Bloch Vector angle
    title_name = "Contrast: " + file_name[13:-8]
    title = "%s"%(title_name)
    label = ['Polar angle $\phi$ [degrees]', 'Avg PMT counts', title]
    def constrast_fit(phi, max_counts, min_counts):
        phi = phi*pi/180
        return (max_counts - min_counts)*sin(phi/2.)**2 + min_counts
    res, res_err = pt.plot_fit(phi, m, constrast_fit, np.array([1000.0, 150.0]),
                               labels=label, axis='default', show=show)

    if res[0] > res[1]:
        max_r = res[0]
        mar_r_err = res_err[0]

        min_r = res[1]
        min_r_err = res_err[1]
    else:
        max_r = res[1]
        mar_r_err = res_err[1]

        min_r = res[0]
        min_r_err = res_err[0]

    contrast = (max_r-min_r)/(max_full-min_full)
    contrast_err = sqrt(res_err[0]**2 + res_err[1]**2) / (max_full-min_full)

    return contrast, contrast_err

def sq_figures(max_c, min_c, N, N_err, sigA, k0, Jbar_1kHz, save=False, extent='on'):

    # Get data
    file_name, scandata, m, pmterr, trials, data = hfGUIdata.get_raw_counts()

    #data = np.genfromtxt(file_name,delimiter=",",names=True,dtype=None)
    arm_time = np.mean(data['arm_t'])
    pi_time = np.mean(data['middle_t'])
    detune = np.mean(data['det_n'])*1e3/(arm_time)  # kHz

    # Calculate derived quantities
    phi = (scandata/pi*180.0) #degrees
    chi = 2*Jbar_1kHz/detune/N  # interaction strength, scaled for detuning
    chi_err = chi*N_err/N

    pmterr_uncert = pmterr/sqrt(2*trials)
    m_avg = np.mean(m)
    K = max_c - min_c

    def noise_model(phi, Ntot, sigA, k0):
        phi = phi*pi/180
        varPN = K**2 * (sin(phi)**2)/4.0/Ntot  # extra factor of K is for conversion to photons
        varAN = (sigA * sin(phi))**2  # sigA is the std dev of rotation of theta in radians
        varSN = K*(sin(phi/2.)**2) + min_c
        varTN = k0**2 * (K*(sin(phi/2.)**2) + min_c)**2
        return sqrt(varPN + varAN + varSN + varTN)

    # store the baselines in tuples to keep uncertainty
    # to predict the SQL, don't include the added noise in the anti-squeezed quad
    NL_sdev = (noise_model(90., N, 0.0, k0), noise_model(90., N+N_err, 0.0, k0), noise_model(90., N-N_err, 0.0, k0))
    PN_sdv = (K/2.0/sqrt(N), K/2.0/sqrt(N+N_err),K/2.0/sqrt(N-N_err))

    ############### Std dev plot ##########################
    title_name = "data set: "+file_name[13:-8]
    title = "%s, $t_{a}$:%d us, $\delta$:%.3g kHz"%(title_name, arm_time, detune)

    label = ['Final rotation [degrees]', 'Std. Dev. PMT counts', title]

    plt.close()
    plt.errorbar(phi, pmterr, yerr=pmterr_uncert, fmt='o')
    plt.fill_between(phi, NL_sdev[0], y2=NL_sdev[1], facecolor=pt.red, alpha=0.5 )
    plt.xlabel(label[0])
    plt.ylabel(label[1])
    plt.title(title)

    # try to calculated the squeezed state from OAT
    psi = np.linspace(0,350,num=200)*pi/180.
    t = 2*arm_time*1e-6

    varJz, delta = ODF.OATjgb(N, Jbar_1kHz, psi, t)

    sig_squeeze = sqrt(varJz/N**2 * K**2)
    sig_squeeze = sqrt(sig_squeeze**2 + 0.0**2 + np.mean(m) + (k0*np.mean(m))**2)
    plt.plot(180/pi*(psi),sig_squeeze,'-',c=pt.blue)
    plt.show()

    ################## dB plot ###########################
    R = pmterr**2/PN_sdv[0]**2
    R_err = 2*pmterr/PN_sdv[0]**2 * pmterr_uncert
    R_err = 10.0*np.log10(R+R_err) - 10*np.log10(R)
    R = 10*np.log10(R)
    SQL = (PN_sdv[1]**2/PN_sdv[0]**2, PN_sdv[2]**2/PN_sdv[0]**2)
    SQL = 10*np.log10(SQL)

    R_pred = 10*np.log10(sig_squeeze**2/PN_sdv[0]**2)
    plt.close()
    ax = plt.axes()
    pt.set_plot_mode(ax,mode='clean')

    plt.plot(180/pi*(psi),R_pred,c=pt.blue)
    plt.errorbar(phi, R, yerr=R_err, fmt='ko')
    plt.fill_between(phi, SQL[0], y2=SQL[1], facecolor=pt.red, alpha=0.5 )
    plt.xlabel(r'Rotation angle $\psi$ [degree]')
    plt.ylabel('Observed spin variance [dB]')
    plt.xticks([0,90,180,270,360])
    plt.yticks([-6,-3,0,3,6,9,12])
    plt.axis(extent)
    #plt.axis(pt.auto_extent(phi,R))
    if save is True:
        plt.savefig(str(int(arm_time))+"_dB_squeeze.png",transparent=True,dpi=300)
    plt.show()

    # save output file of data
    R = pmterr**2/PN_sdv[0]**2
    pts = np.ones(np.size(R))
    int_time = pts*t
    detuning = pts*detune
    N = pts*N
    PSN = m/PN_sdv[0]**2
    output_labels = 'interaction_time_us, detuning_kHz, N, spin_var_rel_PN, tomography_angle_deg, PSN'
    output = np.transpose([int_time, detuning, N, R, phi, PSN])
    #out = np.array([output_labels, output])
    csvname = file_name[13:-8]+"_dataout.csv"
    np.savetxt(csvname, output, delimiter=',', fmt='%s', header=output_labels)

    return 2*arm_time, detune, m, np.min(pmterr), phi, pmterr, psi, sig_squeeze

def data_point_histogram(pn, max_c, min_c, binwidth=0.025, save=False):
    # Get data
    file_name, scandata, counts_data, data = hfGUIdata.get_raw_counts_hist()

    #data = np.genfromtxt(file_name,delimiter=",",names=True,dtype=None)
    arm_time = np.mean(data['arm_t'])
    pi_time = np.mean(data['middle_t'])
    det_n = np.mean(data['det_n'])
    rot_time = np.mean(data['final_t']) * np.ones(np.size(scandata))

    # Calculate derived quantities
    detune = det_n*1e3/(arm_time)  # kHz
    phi = (rot_time/pi_time*180.0)  # degrees

    z_data = (counts_data-min_c)/(max_c - min_c)

    plt.close()
    bin_width = binwidth
    bin_range = np.arange(min(z_data[pn]), max(z_data[pn]) + bin_width, bin_width)
    out = plt.hist(z_data[pn],bins=bin_range)
    #  param = norm.fit(counts_data[pn])
    #  pdf_fitted = dist.pdf(np.arange(1000), *param[:-2], loc=param[-2], scale=param[-1]) * 1000
    x_data = np.array([i*np.ones(np.size(z_data[0])) for i in phi])
    x = x_data.flatten()
    y = z_data.flatten()
    #plt.hist2d(x,y,bins=[25,50])
    #plt.axis([0,180,0,1200])
    plt.xlabel('Counts')
    plt.ylabel('Number of Instances')
    title = 'Angle:{0:.1f} deg, t_a: {1} us, #{2}'.format(phi[pn],arm_time,pn)
    plt.title(title)
    plt.axis([0,1.0,0,np.max(out[0])])

    if save is not False:
        name_out = save+'/phi_'+'{:.0f}'.format(phi[pn])+'_{}'.format(pn)+'.png'
        try:
            plt.savefig(name_out, format='png', bbox='tight')
        except IOError as exception:
            print("No such file or directory")

    return counts_data, x_data

def hist_data_browser(max_c, min_c, num):
    directory = "data_browswer"
    if not os.path.exists(directory):
        os.makedirs(directory)
    for i in range(num):
        data_point_histogram(i, max_c, min_c, save=directory)

def OAT_decoh(psi, tsq, J, N, G_el, G_ud, G_du):
    """
    Model for OAT accounting for decohrence based on paper by MFF
    psi: angle of final rotation about x, radians
    tsq: total interaction tim, seconds
    J: interaction strength from mean field rotation
    N: ion number
    G_el: elastic scattering rate (1/sec)
    G_ud: up down rate (1/sec)
    G_du: down up rate (1/sec)

    Returns:
    Jz_std: normalize to sqrt(N)/2
    C: fringe contrast
    opt_squ_angle: angle of minimum in Jz_std, in radians
    """

    G_r = G_ud + G_du
    G = (G_el + G_r)/2.

    # Key functions
    def Phi(t, J, G_ud, G_du):
        theta = pi/2.
        g = (G_ud-G_du)/4.
        G_r = G_ud + G_du
        s = 2*J/N + 2j*g  # this is the param in the paper where N should go
        r= G_ud*G_du

        a = t*sqrt(s**2-r)

        out = exp(-G_r*t/2.)*( cos(a) + 0.5*(G_r+4j*J*cos(theta))*t*np.sinc(a) )
        return out

    def Psi(t, J, G_ud, G_du):
        theta = pi/2.
        g = (G_ud-G_du)/4.
        G_r = G_ud + G_du
        s = 2*J/N + 2j*g
        r= G_ud*G_du

        a = t*sqrt(s**2-r)

        out = exp(-G_r*t/2.)*( cos(theta)*cos(a) + 0.5*(2j*s - 4*g - G_r*cos(theta))*t*np.sinc(a) )
        return out

    # collective spin length
    Sx = 0.5*exp(-G*tsq) * (Phi(tsq,J,G_ud, G_du))**(N-1)
    Sx = np.real(Sx)

    # pm correlations
    Cmz = 0.5*exp(-G*tsq) * Psi(tsq, -J ,G_ud, G_du) * (Phi(tsq, -J, G_ud, G_du))**(N-2)
    Cpz = 0.5*exp(-G*tsq) * Psi(tsq, J, G_ud, G_du) * (Phi(tsq, J, G_ud, G_du))**(N-2)

    Cpp = 0.25*exp(-2*G*tsq) * (Phi(tsq, 2*J, G_ud, G_du))**(N-2)
    Cmm = 0.25*exp(-2*G*tsq) * (Phi(tsq, -2*J, G_ud, G_du))**(N-2)
    Cpm = 0.25*exp(-2*G*tsq) * (Phi(tsq, 0.0, G_ud, G_du))**(N-2)
    Cmp = Cpm

    # coorediate correlations
    CyyA = -0.25*(Cpp+Cmm-Cpm-Cmp)
    CyzA = -(0.25j)*(Cpz-Cmz)
    CzyA = CyzA

    SyyA = N*(N-1)*CyyA + N/4.
    SyzA = N*(N-1)*CyzA + 0.5j*N*Sx
    SzyA = N*(N-1)*CzyA - 0.5j*N*Sx
    SzzA = N/4.

    # calc squeezing
    z = (SzyA+SyzA)/(SzzA-SyyA)
    opt_squ_angle = np.real(0.5*np.arctan(z))
    Jz_std = sqrt(cos(-psi)**2*SzzA + sin(-psi)**2*SyyA + sin(-psi)*cos(-psi)*(SyzA+SzyA))
    R = (Jz_std/(sqrt(N)/(2.0)))**2  # reduction in spin noise variance
    C = 2*Sx  # fringe contrast

    Xi2 = R * (C)**(-2)  # Ramsey squeezing parameter squared

    return np.array([Jz_std, C, opt_squ_angle])

