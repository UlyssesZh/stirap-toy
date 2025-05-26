#!/usr/bin/env python

from numpy import linspace, loadtxt, array
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

plot_dir = 'paper/amplitude-sweep'
output_file = 'paper/amplitude-sweep.csv'
output_plot = 'paper/amplitude-sweep.pdf'
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
	alpha0_vals = []
	n1_vals = []
	n2_vals = []
else:
	alpha0_vals = list(data[:, 1])
	n1_vals = list(data[:, 2])
	n2_vals = list(data[:, 3])

for i1, alpha0 in enumerate(linspace(500, 3500, 51)):
	if i1 < len(alpha0_vals):
		continue

	param.alpha0 = alpha0
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
		plt.title(f'$\\alpha_0={alpha0:.3f}$')
		plt.savefig(f'{plot_dir}/{i1:03d}.pdf', transparent=True, format='pdf', bbox_inches='tight')
		plt.close()
		
	alpha0_vals.append(alpha0)
	n1_val = output.expect[1][-1]
	n2_val = output.expect[2][-1]
	n1_vals.append(n1_val)
	n2_vals.append(n2_val)
	with open(output_file, 'a') as f:
		f.write(f'{i1},{alpha0},{n1_val},{n2_val}\n')

plt.figure()
plt.plot(alpha0_vals, n2_vals)
plt.xlabel(r"$\alpha_0$")
plt.ylabel(r"$\left<n_2\right>$")
plt.savefig(output_plot, transparent=True, format='pdf', bbox_inches='tight')
