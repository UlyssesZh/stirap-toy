#qutip imports
from qutip import *
from qutip.ui.progressbar import BaseProgressBar, TextProgressBar
from qutip.tensor import flatten

#standard python imports
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.constants as sc

#Our imports
from pumps import *

#phi are the total phase rotations of the a*b.dag() and a.dag()*b operators
def Hamiltonian(params_dict, fig, ax1):

    #assign necessary paremeters from the dictionary
    N, Nopt, omega_1, omega_2, omega_c, g1, g2,alpha1, alpha2,  sigma_1, sigma_2, gamma_1, gamma_2, phi_1, phi_2, phi_1_p, phi_2_p, phi_spring_1_plus,phi_spring_1_minus,phi_spring_2_minus,phi_spring_2_plus, delta, tau, t_delay, t_min, t_max, numstep, tlist, n = params_dict["N"], params_dict["Nopt"], params_dict["omega_1"], params_dict["omega_2"], params_dict["omega_c"], params_dict["g1"], params_dict["g2"], params_dict["alpha1"], params_dict["alpha2"], params_dict["sigma_1"], params_dict["sigma_2"],params_dict["gamma_1"],params_dict["gamma_2"],params_dict["phi_1"], params_dict["phi_2"], params_dict["phi_1_p"], params_dict["phi_2_p"], params_dict["phi_spring_1+"],params_dict["phi_spring_1-"], params_dict["phi_spring_2-"], params_dict["phi_spring_2+"], params_dict["delta"], params_dict["tau"], params_dict["t_delay"], params_dict["t_min"], params_dict["t_max"], params_dict["num_step"], params_dict["tlist"], params_dict["n"]
    #desctruction operators on the tensor product space
    a  = tensor(destroy(Nopt), qeye(N), qeye(N))
    b1 = tensor(qeye(Nopt), destroy(N), qeye(N))
    b2 = tensor(qeye(Nopt), qeye(N), destroy(N))
    d  = tensor(qeye(Nopt), qeye(N), qeye(N))

    #position and momentum operators via ladder operators
    # x1 = x1zpf * (b1 + b1.dag())
    # x2 = x2zpf * (b2 + b2.dag())
    # p1 = -1j*meff1*omega_1*x1zpf*(b1 - b1.dag())
    # p1 = -1j*meff2*omega_2*x2zpf*(b2 - b2.dag())

    #the Hamiltonian
    H11a_1 = sc.hbar * g1 * alpha1 * a * (b1.dag())
    H11a_2 = sc.hbar * g1 * alpha1 * a * (b1)

    H11b_1 = sc.hbar * g1 * alpha1 * a.dag() * (b1.dag())
    H11b_2 = sc.hbar * g1 * alpha1 * a.dag() * (b1)

    H12a_1 = sc.hbar * g1 * alpha2 * a * (b1.dag())
    H12a_2 = sc.hbar * g1 * alpha2 * a * (b1)


    H12b_1 = sc.hbar * g1 * alpha2 * a.dag() * (b1.dag())
    H12b_2 = sc.hbar * g1 * alpha2 * a.dag() * (b1)


    H21a_1 = sc.hbar * g2 * alpha1 * a * (b2.dag())
    H21a_2 = sc.hbar * g2 * alpha1 * a * (b2)

    H21b_1 = sc.hbar * g2 * alpha1 * a.dag() * (b2.dag())
    H21b_2 = sc.hbar * g2 * alpha1 * a.dag() * (b2)

    H22a_1 = sc.hbar * g2 * alpha2 * a * (b2.dag())
    H22a_2 = sc.hbar * g2 * alpha2 * a * (b2)

    H22b_1 = sc.hbar * g2 * alpha2 * a.dag() * (b2.dag())
    H22b_2 = sc.hbar * g2 * alpha2 * a.dag() * (b2)

        #H12a = sc.hbar * g1 * alpha2 * a * (b1.dag())
        #H12b = sc.hbar * g1 * alpha2 * a.dag() * (b1)
        #H12a = a - a
        #H12b = a - a
        #H21a = sc.hbar * g2 * alpha1 * a * (b2.dag())
        #H21b = sc.hbar * g2 * alpha1 * a.dag() * (b2)
        #H21a = a - a
        #H21b = a - a


    Hm = omega_1 * b1.dag() * b1 + omega_2 * b2.dag()*b2
    Hm = Hm - Hm
    #The pulses you would like to implement
    #pulse_1() = fSTIRAP_PULSE(tlist,  tau, t_delay, sigma_1) *  reversed_fSTIRAP_PULSE(tlist,  tau, t_delay, sigma_1)

    #The total hamiltonian, that encodes the time dependence by multiplying the cavity annihalation/creation operators by a plane wave
#    H = [
#    Hm,
##    [H11a_1,phase(tlist,-1, phi_1) * gaussian1(tlist, tau, t_delay, sigma_1)],
#    [H11a_2,phase(tlist,-1, phi_1) * gaussian1(tlist, tau, t_delay, sigma_1)],
#    [H11b_1,phase(tlist,+1, phi_1_p) * gaussian1(tlist, tau, t_delay, sigma_1)],
##    [H11b_2,phase(tlist,+1, phi_1) * gaussian1(tlist, tau, t_delay, sigma_1)],
#    [H12a_1,phase(tlist,-1, phi_1_p) * gaussian2(tlist, tau, t_delay, sigma_1)],
#    [H12a_2,phase(tlist,-1, phi_1) * gaussian2(tlist, tau, t_delay, sigma_1)],
#    [H12b_1,phase(tlist,+1, phi_1_p) * gaussian2(tlist, tau, t_delay, sigma_1)],
#    [H12b_2,phase(tlist,+1, phi_1) * gaussian2(tlist, tau, t_delay, sigma_1)],
#    [H21a_1,phase(tlist,-1, phi_2_p) * gaussian1(tlist, tau, t_delay, sigma_1)],
#    [H21a_2,phase(tlist,-1, phi_2) * gaussian1(tlist, tau, t_delay, sigma_1)],
#    [H21b_1,phase(tlist,+1, phi_2_p) * gaussian1(tlist, tau, t_delay, sigma_1)],
#    [H21b_2,phase(tlist,+1, phi_2) * gaussian1(tlist, tau, t_delay, sigma_1)],
##    [H22a_1,phase(tlist,-1, phi_2) * gaussian2(tlist, tau, t_delay, sigma_1)],
#    [H22a_2,phase(tlist,-1, phi_2) * gaussian2(tlist, tau, t_delay, sigma_1)],
#    [H22b_1,phase(tlist,+1, phi_2_p) * gaussian2(tlist, tau, t_delay, sigma_1)],
##    [H22b_2,phase(tlist,+1, phi_2) * gaussian2(tlist, tau, t_delay, sigma_1)]
#    ]
#last four was ++--

#    H = [
#    Hm,
#    [H11a_1,phase(tlist,-1,0) * gaussian1(tlist, tau, t_delay, sigma_1)],
#    [H11a_2,(phase(tlist,-1,2*omega_1)) * gaussian1(tlist, tau, t_delay, sigma_1)],
#    [H11b_1,(phase(tlist,+1,2*omega_1)) * gaussian1(tlist, tau, t_delay, sigma_1)],
#    [H11b_2,phase(tlist,+1,0) * gaussian1(tlist, tau, t_delay, sigma_1)],
#    [H12a_1,phase(tlist,-1,omega_1-omega_2) * gaussian2(tlist, tau, t_delay, sigma_1)],
#    [H12a_2,phase(tlist,-1,omega_2+omega_1) * gaussian2(tlist, tau, t_delay, sigma_1)],
#    [H12b_1,phase(tlist,+1,omega_1+omega_2) * gaussian2(tlist, tau, t_delay, sigma_1)],
#    [H12b_2,phase(tlist,+1,omega_2-omega_1) * gaussian2(tlist, tau, t_delay, sigma_1)],
#    [H21a_1,phase(tlist,-1,omega_1-omega_2) * gaussian1(tlist, tau, t_delay, sigma_1)],
#    [H21a_2,phase(tlist,-1,omega_1+omega_2) * gaussian1(tlist, tau, t_delay, sigma_1)],
#    [H21b_1,phase(tlist,+1,omega_1+omega_2) * gaussian1(tlist, tau, t_delay, sigma_1)],
#    [H21b_2,phase(tlist,+1,omega_1-omega_2) * gaussian1(tlist, tau, t_delay, sigma_1)],
#    [H22a_1,phase(tlist,-1,0) * gaussian2(tlist, tau, t_delay, sigma_1)],
#    [H22a_2,(phase(tlist,-1,2*omega_2)) * gaussian2(tlist, tau, t_delay, sigma_1)],
#    [H22b_1,(phase(tlist,+1,2*omega_2)) * gaussian2(tlist, tau, t_delay, sigma_1)],
#    [H22b_2,phase(tlist,+1,0) * gaussian2(tlist, tau, t_delay, sigma_1)]
#    ]

    H = [
    Hm,
    [H11a_1,phase(tlist,-1,omega_1 - omega_1) * gaussian1(tlist, tau, t_delay, sigma_1)],
    [H11a_2,(phase(tlist,-1,omega_1+omega_1)) * gaussian1(tlist, tau, t_delay, sigma_1)],
    [H11b_1,(phase(tlist,+1,omega_1+omega_1)) * gaussian1(tlist, tau, t_delay, sigma_1)],
    [H11b_2,phase(tlist,+1,omega_1 - omega_1) * gaussian1(tlist, tau, t_delay, sigma_1)],
    [H12a_1,phase(tlist,-1,n*omega_2-omega_1) * gaussian2(tlist, tau, t_delay, sigma_1)],
    [H12a_2,phase(tlist,-1,n*omega_2+omega_1) * gaussian2(tlist, tau, t_delay, sigma_1)],
    [H12b_1,phase(tlist,+1,n*omega_2+omega_1) * gaussian2(tlist, tau, t_delay, sigma_1)],
    [H12b_2,phase(tlist,+1,n*omega_2-omega_1) * gaussian2(tlist, tau, t_delay, sigma_1)],
    [H21a_1,phase(tlist,-1,omega_1-omega_2) * gaussian1(tlist, tau, t_delay, sigma_1)],
    [H21a_2,phase(tlist,-1,omega_1+omega_2) * gaussian1(tlist, tau, t_delay, sigma_1)],
    [H21b_1,phase(tlist,+1,omega_1+omega_2) * gaussian1(tlist, tau, t_delay, sigma_1)],
    [H21b_2,phase(tlist,+1,omega_1-omega_2) * gaussian1(tlist, tau, t_delay, sigma_1)],
    [H22a_1,phase(tlist,-1,omega_2 - n*omega_2) * gaussian2(tlist, tau, t_delay, sigma_1)],
    [H22a_2,(phase(tlist,-1,omega_2 + n*omega_2)) * gaussian2(tlist, tau, t_delay, sigma_1)],
    [H22b_1,(phase(tlist,+1,omega_2 + n*omega_2)) * gaussian2(tlist, tau, t_delay, sigma_1)],
    [H22b_2,phase(tlist,+1,omega_2 - n*omega_2) * gaussian2(tlist, tau, t_delay, sigma_1)]
    ]


    #Plot the pumps
    ax2 = ax1[0].twinx()
    ax3 = ax1[0].twinx()

    axes = [ax1[0], ax2, ax3]

    fig.subplots_adjust(right=0.75)
    ax3.spines['right'].set_position(('axes', 1.2))
    ax3.set_frame_on(True)
    ax3.patch.set_visible(False)

    #ax1[0].legend(loc = 'upper right')
    ax3.legend(loc = 'upper center')
    #ax1[1].legend(loc = 'upper right')

    ax2.plot(tlist, gaussian1(tlist, tau, t_delay, sigma_1), marker = '', color = 'red',label = 'Constant Pulse')
    ax2.plot(tlist, gaussian2(tlist, tau, t_delay, sigma_1), marker = '', color = 'blue', label = 'P Pulse')
#    ax3.plot(tlist,g1 * alpha1 * gaussian(tlist, tau, t_delay, 0.05), marker = '', color = 'orange', label = 'Spring pulse')
#    ax3.tick_params(axis='y', colors = 'orange')
    return [H, a, b1, b2, d, fig, ax1]
