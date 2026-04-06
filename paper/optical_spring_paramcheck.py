

from qutip import *

from qutip.ui.progressbar import BaseProgressBar, TextProgressBar
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.constants as sc

from qutip.tensor import flatten
#from IPython.display import Image




def get_gamma_opt(alpha, g, kappa, delta, omega):
    #optical spring effect decay rate
    gamma_opt = (alpha* g)**2 * (kappa / ((kappa**2) / 4 + (delta + omega)**2) - (kappa / ((kappa**2)/4+(delta - omega)**2)))
    #gamma_opt = g**2 * omega * (-8 * omega * kappa) / (kappa**2 / 4 + 4*omega)**2

    print("gamma opt = {}".format(gamma_opt))
    return gamma_opt
    #optical spring effect frequency shift

def get_delta_omega(alpha, g, kappa, delta, omega):
    delta_omega = (alpha * g)**2*((delta - omega) / ((kappa**2) /4 + (delta - omega)**2) + ((delta + omega) / ((kappa**2) /4 + (delta + omega)**2)))
    print("delta omega = {}".format(delta_omega))
    return delta_omega


def plot_gamma_opt_decay(tlist, gamma_opt, t_min):
    profile = np.exp(-gamma_opt * (tlist - t_min))
    return profile

def plot_real_spring_phase(tlist, delta_omega, t_min):
    profile = np.real(np.exp(-1j *0.5* delta_omega * tlist))
    return profile

def get_gamma_opt1(alpha, g, kappa, delta, omega):
    #optical spring effect decay rate
    gamma_opt1 = (alpha* g)**2 * (kappa / ((kappa**2) / 4 + (delta + omega)**2) - (kappa / ((kappa**2)/4+(delta - omega)**2)))
    #gamma_opt = g**2 * omega * (-8 * omega * kappa) / (kappa**2 / 4 + 4*omega)**2

    print("gamma opt 1 = {}".format(gamma_opt1))
    return gamma_opt1


def get_gamma_opt2(alpha, g, kappa, delta, omega):
    #optical spring effect decay rate
    gamma_opt2 = (alpha* g)**2 * (kappa / ((kappa**2) / 4 + (delta + omega)**2) - (kappa / ((kappa**2)/4+(delta - omega)**2)))
    #gamma_opt = g**2 * omega * (-8 * omega * kappa) / (kappa**2 / 4 + 4*omega)**2

    print("gamma opt 2 = {}".format(gamma_opt2))
    return gamma_opt2

def get_delta_omega1(alpha, g, kappa, delta, omega):
    delta_omega1 = (alpha * g)**2*((delta - omega) / ((kappa**2) /4 + (delta - omega)**2) + ((delta + omega) / ((kappa**2) /4 + (delta + omega)**2)))
    print("delta omega 1 = {}".format(delta_omega1))
    return delta_omega1

def get_delta_omega2(alpha, g, kappa, delta, omega):
    delta_omega2 = (alpha * g)**2*((delta - omega) / ((kappa**2) /4 + (delta - omega)**2) + ((delta + omega) / ((kappa**2) /4 + (delta + omega)**2)))
    print("delta omega 2 = {}".format(delta_omega2))
    return delta_omega2

def get_gamma_optical_cooling(g, alpha, kappa):
    gamma_optical_cooling = 4 * (g * alpha)**2 / (kappa)
    print("gamma_optical_cooling = {}".format(gamma_optical_cooling))
    return gamma_optical_cooling

def mean_phonon_numer_optical_cooling(tlist, gamma_optical_cooling, t_min, nth, kappa, gamma, omega_m):
    mean_phonon_numer_optical_cooling = nth * (gamma + gamma_optical_cooling * np.exp(-gamma_optical_cooling * (tlist-t_min)))/(gamma + gamma_optical_cooling) + (kappa**2/(16*omega_m**2)) * (1-np.exp(-gamma_optical_cooling*(tlist-t_min)))
    return mean_phonon_numer_optical_cooling
