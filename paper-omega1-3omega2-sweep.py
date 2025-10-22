#!/usr/bin/env python

from numpy import linspace, pi, exp, sqrt, sin, cos, loadtxt, array, argmax
import matplotlib.pyplot as plt
from qutip import mesolve, fidelity, hinton, plot_wigner

import dim
import param
import t_dep
from op import setup

import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

plot_dir = 'paper/omega1-3omega2-sweep'
output_file = 'paper/omega1-3omega2-sweep.csv'
output_plot = 'paper/omega1-3omega2-sweep.pdf'
import os
if not os.path.exists(plot_dir):
	os.makedirs(plot_dir)

if os.path.exists(output_file):
	data = loadtxt(output_file, delimiter=',')
else:
	data = []
	with open(output_file, 'w') as f:
		pass
if len(data) == 0:
	omega1_vals = []
	omega2_vals = []
	n1_vals = []
	n2_vals = []
else:
	omega1_vals = list(data[:, 1])
	omega2_vals = list(data[:, 2])
	n1_vals = list(data[:, 3])
	n2_vals = list(data[:, 4])

plt.rcParams.update({
	'text.latex.preamble': r'\usepackage{lmodern}\usepackage{siunitx}',
	'text.usetex': True,
	'font.size': 11,
	'font.family': 'lmodern',
	'lines.linewidth': 1
})

count = 0
for i, omega1 in enumerate(linspace(0.02e6-1800, 0.02e6+1800, 401)*2*pi):
	omega2 = (0.06e6 - (omega1/(2*pi)-0.02e6)/3)*2*pi
	count += 1
	if count <= len(data):
		continue

	param.omega1 = omega1
	param.omega2 = omega2
	H, c_ops, e_ops, psi0, psi1 = setup()

	options={
		#'nsteps': 2000000000,
		'progress_bar': 'tqdm',
		'store_final_state': True
	}

	t, alpha_n1, alpha_n2 = t_dep.t, t_dep.alpha_n1, t_dep.alpha_n2
	output = mesolve(H, psi0, t, c_ops, e_ops=e_ops, options=options)
	final = output.final_state
	fid = fidelity(psi1, final)
	print(f'Fidelity: {fid**2}')
	
	plt.figure()
	plt.plot(t*1000, output.expect[1], label='Mechanical 1 num expec')
	plt.plot(t*1000, output.expect[2], label='Mechanical 2 num expec')
	plt.plot(t*1000, alpha_n1, label="Pulse 1")
	plt.plot(t*1000, alpha_n2, label="Pulse 2")
	plt.legend()
	plt.xlabel("Time (ms)")
	plt.title(r'$\omega_1/2\pi=\SI{%.6f}{MHz},\omega_2/2\pi=\SI{%.6f}{MHz}$' % (omega1*1e-6/(2*pi), omega2*1e-6/(2*pi)))
	plt.savefig(f'{plot_dir}/{i}.pdf', transparent=True, format='pdf', bbox_inches='tight')
	plt.close()

	omega1_vals.append(omega1)
	omega2_vals.append(omega2)
	n1_val = output.expect[1][-1]
	n2_val = output.expect[2][-1]
	n1_vals.append(n1_val)
	n2_vals.append(n2_val)
	with open(output_file, 'a') as f:
		f.write(f'{i},{omega1},{omega2},{n1_val},{n2_val}\n')

plt.xlabel(r'$\omega_1/2\pi$ (MHz)')
plt.ylabel(r'$\left<n_2\right>$')
plt.title(r'$\omega_2/2\pi-\SI{0.06}{MHz}=-\left(\omega_1/2\pi-\SI{0.02}{MHz}\right)/3$')
plt.plot(array(omega1_vals)*1e-6/(2*pi), n2_vals)
plt.savefig(output_plot, transparent=True, format='pdf', bbox_inches='tight')
