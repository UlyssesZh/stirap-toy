#qutip imports
from qutip import *
from qutip.ui.progressbar import BaseProgressBar, TextProgressBar
from qutip.tensor import flatten

#standard python imports
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.constants as sc


#plane wave
def phase(tlist,  sign, phi):
    return np.exp(sign *1j *tlist * phi)

#gaussian profiles for the pumps, with two terms in the s pump for partial stirap
def gaussian(tlist,  tau, t_delay, sigma):
    #profile = np.sqrt(2) / 2 * np.exp(-1 * pow((((tlist + t_delay) - tau)/sigma_1),2)) + np.sqrt(2) / 2 * np.exp(-1 * pow(((-(tlist - t_delay) - tau)/sigma_1),2))
    profile = 1/2.2 * np.exp(-1 * pow((((tlist))/sigma),2))
    #profile =np.sqrt(2) / 2 * np.exp(-1 * pow((((tlist + t_delay) - tau)/sigma),2)) + np.sqrt(2) / 2 * np.exp(-1 * pow(((-(tlist - t_delay) - tau)/sigma),2))
    return profile

def constant(tlist, delta):
    if -1e-1 <= tlist <= 0.0:
        return tlist - tlist + delta
    elif 0.0<= tlist <= 1e-1:
        return tlist - tlist - delta
    return 0.0
constant_vector = np.vectorize(constant)
def constant_blue(tlist, g1, alpha1):
    if -1e-1 <= tlist <= 0.0:
        return tlist - tlist + 1
    elif 0.0<= tlist <= 0.5e-1:
        return tlist - tlist + 1
    return 0.0
constant_vector_blue = np.vectorize(constant_blue)

def fSTIRAP_PULSE(tlist,  tau, t_delay, sigma):
    profile =np.sqrt(2) / 2 * np.exp(-1 * pow((((tlist - t_delay) - tau)/sigma),2))
    return profile

def fSTIRAP_P(tlist,  tau, t_delay, sigma):
    profile =(np.sqrt(2) / 2) * np.exp(-1 * pow((((tlist + t_delay) - tau)/sigma),2))
    return profile

def fSTIRAP_S(tlist,  tau, t_delay, sigma):
    profile = np.exp(-1 * pow((((tlist+t_delay) + tau)/sigma),2)) + (np.sqrt(2) / 2) * np.exp(-1 * pow((((tlist+t_delay) - tau)/sigma),2))
    return profile

def ffSTIRAP_S(tlist,  tau, t_delay, sigma):
    profile = np.exp(-1 * pow((((tlist+t_delay) - tau)/sigma),2))
    return profile

def reversed_fSTIRAP_P(tlist,  tau, t_delay, sigma):
    profile =np.sqrt(2) / 2 * np.exp(-1 * pow((((-1*tlist + t_delay) - tau)/sigma),2))
    return profile

def reversed_fSTIRAP_PULSE(tlist,  tau, t_delay, sigma):
    profile =np.sqrt(2) / 2 * np.exp(-1 * pow(((-(tlist - t_delay) - tau)/sigma),2))
    return profile

def gaussian1(tlist,  tau, t_delay, sigma):
    profile =np.exp(-1 * pow((((tlist-tau+t_delay))/sigma),2))
    return profile

def gaussian2(tlist,  tau, t_delay, sigma):
    profile = np.exp(-1 * pow((((tlist+tau))/sigma),2))
    return profile

def half_pi_pulse(tlist, tau_half_pi, t_delay, sigma):
    profile = np.exp(-1 * pow((((tlist))/tau_half_pi),2)/2)
    return profile

def constant_pulse(tlist):
    profile = tlist - tlist + 1
    return profile

def constant_red(tlist, delta):
    if -1e-1 <= tlist <= 1e-1:
        return delta
    elif 1e-1 <= tlist <= 4e-1:
        return 0.0
    return 0.0

def constant_blue(tlist, delta):
    if -1e-1 <= tlist <= 1e-1:
        return 0.0
    elif 1e-1 <= tlist <= 4e-1:
        return -1 * delta
    return 0.0

constant_blue_vector = np.vectorize(constant_blue)
constant_red_vector = np.vectorize(constant_red)
