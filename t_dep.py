from numpy import linspace, exp, sin, cos, sign

import param

import sys
this = sys.modules[__name__]

num_sigma = 5.0
nstep_per_sec = 3.75e6

def setup():
	param.setup()
	tau1, tau2, sigma1, sigma2 = param.tau1, param.tau2, param.sigma1, param.sigma2
	theta, alpha0 = param.theta, param.alpha0

	boundaries = [
		tau1 + num_sigma*sigma1,
		tau1 - num_sigma*sigma1,
		tau2 + num_sigma*sigma2,
		tau2 - num_sigma*sigma2
	]
	this.t_min = min(boundaries)
	this.t_max = max(boundaries)
	this.t = linspace(t_min, t_max, int(nstep_per_sec*(t_max-t_min)))
	this.pulse1 = exp(-(t-tau1)**2/sigma1**2)
	this.pulse2 = exp(-(t-tau2)**2/sigma2**2)
	this.alpha_n1 = sin(theta)*pulse1
	this.alpha_n2 = pulse2 + cos(theta)*pulse1
	this.alpha1 = alpha0*alpha_n1
	this.alpha2 = alpha0*alpha_n2
