from numpy import linspace, exp, sin, cos

from param import *

t_min = -2e-3
t_max = 2e-3
nstep = 15000

t = linspace(t_min, t_max, nstep)

pulse1 = exp(-(t-tau1)**2/sigma1**2)
pulse2 = exp(-(t-tau2)**2/sigma2**2)
alpha_n1 = sin(theta)*pulse1
alpha_n2 = pulse2 + cos(theta)*pulse1
alpha1 = alpha0*alpha_n1
alpha2 = alpha0*alpha_n2
