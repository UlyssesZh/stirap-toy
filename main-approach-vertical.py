#!/usr/bin/env python

from numpy import linspace, pi, exp, sqrt, sin, cos, concatenate
import matplotlib.pyplot as plt
from qutip import mesolve, fidelity, hinton, plot_wigner

import dim
import param
import t_dep
from op import setup

import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)


plt.rcParams.update({
	'text.latex.preamble': r'\usepackage{lmodern}\usepackage{siunitx}',
	'text.usetex': True,
	'font.size': 11,
	'font.family': 'lmodern',
	'lines.linewidth': 1
})

values = concatenate((
	linspace(1.845e6, 1.858e6, 51, endpoint=False),
	linspace(1.858e6, 1.87e6, 11, endpoint=False),
	linspace(1.87e6, 1.88e6, 51, endpoint=False),
	linspace(1.88e6, 1.899e6, 11, endpoint=False),
	linspace(1.899e6, 1.905e6, 21)
))
for i1, omega1 in enumerate(values*2*pi):
	param.omega1 = omega1
	H, c_ops, e_ops, psi0, psi1 = setup()

	options={
		#'nsteps': 2000000000,
		'progress_bar': 'tqdm',
		'store_final_state': True
	}

	t, alpha_n1, alpha_n2 = t_dep.t, t_dep.alpha_n1, t_dep.alpha_n2
	output = mesolve(H, psi0, t, c_ops, e_ops, options=options)
	final = output.final_state
	fid = fidelity(psi1, final)
	print(f'Fidelity: {fid**2}')
	
	if i1%1 == 0:
		plt.figure()
		plt.plot(t*1000, output.expect[1], label='Mechanical 1 num expec')
		plt.plot(t*1000, output.expect[2], label='Mechanical 2 num expec')
		plt.plot(t*1000, alpha_n1, label="Pulse 1")
		plt.plot(t*1000, alpha_n2, label="Pulse 2")
		plt.ylim(0, 2)
		plt.xlabel("Time (ms)")
		plt.ylabel("Amplitude")
		plt.title(f'$\\omega_1/2\pi=\\SI{{{omega1/(2*pi)/1e6:.3f}}}{{MHz}}$')
		plt.savefig(f'plot-approach-vertical/{i1:03d}.pdf', transparent=True, format='pdf', bbox_inches='tight')

	with open('output-approach-vertical.csv', 'a') as f:
		f.write(f'{omega1},{fid}\n')
