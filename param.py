from numpy import pi

from utils import bose

import sys
this = sys.modules[__name__]

omega1 = 1.2e6 * 2*pi
omega2 = 1.8e6 * 2*pi
kappa = 50e3/25 * 2*pi
Delta1 = omega1
Delta2 = omega2
sigma1 = 6e-4
sigma2 = 6e-4
tau1 = sigma1/1.43
tau2 = -tau1
alpha0 = 2000
g1 = 2.5 * 2*pi
g2 = g1
theta = pi/2
T1 = 1e-2
T2 = T1
nth1 = bose(omega1, T1)
nth2 = bose(omega2, T2)
Q1 = 1e9
Q2 = Q1
Gamma_m1 = omega1/Q1
Gamma_m2 = omega2/Q2

def setup():
	this.Delta1 = omega1
	this.Delta2 = omega2
	this.tau2 = -tau1
	this.T2 = T1
	this.nth1 = bose(omega1, T1)
	this.nth2 = bose(omega2, T2)
	this.Q2 = Q1
	this.Gamma_m1 = omega1/Q1
	this.Gamma_m2 = omega2/Q2
