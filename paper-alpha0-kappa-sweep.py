#!/usr/bin/env python

from numpy import linspace, loadtxt, array, pi, argmax
import matplotlib.pyplot as plt
from qutip import mesolve, fidelity, hinton, plot_wigner

import qutip

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

plot_dir = 'paper/alpha0-kappa-sweep'
output_file = 'paper/alpha0-kappa-sweep.csv'
output_plot = 'paper/alpha0-kappa-sweep.pdf'
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
	kappa_vals = []
	n1_vals = []
	n2_vals = []
else:
	alpha0_vals = list(data[:, 2])
	kappa_vals = list(data[:, 3])
	n1_vals = list(data[:, 4])
	n2_vals = list(data[:, 5])

count = 0
for i1, alpha0 in enumerate(linspace(500, 3500, 51)):
	for i2, kappa in enumerate(linspace(0, 10e3*2*pi, 51)):
		count += 1
		if count <= len(data):
			continue

		param.alpha0 = alpha0
		param.kappa = kappa
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
			plt.plot(t*1000, output.expect[0], label='Mechanical 1 prob')
			plt.plot(t*1000, output.expect[1], label='Mechanical 2 prob')
			plt.plot(t*1000, alpha_n1, label="Pulse 1")
			plt.plot(t*1000, alpha_n2, label="Pulse 2")
			plt.ylim(0, 2)
			plt.xlabel("Time (ms)")
			plt.title(r'$\alpha_0=%.3f,\kappa/2\pi=\SI{%.3f}{kHz}$' % (alpha0, kappa/(2*pi)))
			plt.savefig(f'{plot_dir}/{i1:03d}-{i2:03d}.pdf', transparent=True, format='pdf', bbox_inches='tight')
			plt.close()

		alpha0_vals.append(alpha0)
		kappa_vals.append(kappa)
		n1_val = output.expect[1][-1]
		n2_val = output.expect[2][-1]
		n1_vals.append(n1_val)
		n2_vals.append(n2_val)
		with open(output_file, 'a') as f:
			f.write(f'{i1},{i2},{alpha0},{kappa},{n1_val},{n2_val}\n')

x = array(alpha0_vals).reshape(51,51)
y = array(kappa_vals).reshape(51,51) *1e-3/(2*pi)
z = array(n2_vals).reshape(51,51)

max_i = argmax(z)
max_x = x.flat[max_i]
max_y = y.flat[max_i]
print(r'Max <n2>: %.3f at %.3f MHz and %.3f MHz' % (z.flat[max_i], max_x, max_y))

fig, ax = plt.subplots()
mesh = ax.pcolormesh(x, y, z)
mesh.set_edgecolor('face') # https://stackoverflow.com/a/27096694
ax.contour(x, y, z, levels=[0.80, 0.95, 0.99], cmap='copper', linewidths=1)
#plt.plot(max_x, max_y, 'ro')

ax.set_xlabel(r'$\alpha_0$')
ax.set_ylabel(r'$\kappa/2\pi$ (kHz)')

fig.colorbar(mesh, ax=ax).set_label(r'$\left<n_2\right>$')
plt.savefig(output_plot, transparent=True, format='pdf', bbox_inches='tight')
