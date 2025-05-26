#!/usr/bin/env python

from numpy import linspace, pi, exp, sqrt, sin, cos, zeros, save
import matplotlib.pyplot as plt
from qutip import mesolve, fidelity, hinton, plot_wigner

import dim
import param
import t_dep
from op import setup

import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

result = zeros((201, 201, dim.n1, dim.n2))
for i1, omega1 in enumerate(linspace(1.875e6-0.02e6, 1.875e6+0.02e6, 201)*2*pi):
	for i2, omega2 in enumerate(linspace(1.875e6-0.02e6, 1.875e6+0.02e6, 201)*2*pi):
		param.omega1 = param.Delta1 = omega1
		param.omega2 = param.Delta2 = omega2
		H, c_ops, e_ops, psi0, psi1 = setup()

		options={
			#'nsteps': 2000000000,
			'progress_bar': 'tqdm',
			'store_final_state': True
		}
		e_ops = []

		t, alpha_n1, alpha_n2 = t_dep.t, t_dep.alpha_n1, t_dep.alpha_n2
		output = mesolve(H, psi0, t, c_ops, e_ops, options=options)
		final = output.final_state
		result[i1, i2] = final.ptrace((1,2)).diag().reshape((dim.n1, dim.n2))
	save('main-zoomed.npy', result)
