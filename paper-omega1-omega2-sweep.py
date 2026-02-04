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

plot_dir = 'paper/omega1-omega2-sweep'
output_file = 'paper/omega1-omega2-sweep.csv'
output_plot = 'paper/omega1-omega2-sweep.pdf'
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
	omega1_vals = list(data[:, 2])
	omega2_vals = list(data[:, 3])
	n1_vals = list(data[:, 4])
	n2_vals = list(data[:, 5])

plt.rcParams.update({
	'text.latex.preamble': r'\usepackage{lmodern}\usepackage{siunitx}',
	'text.usetex': True,
	'font.size': 11,
	'font.family': 'lmodern',
	'lines.linewidth': 1
})

count = 0
for i1, omega1 in enumerate(linspace(0, 1.5e6, 151)[1:]*2*pi):
	for i2, omega2 in enumerate(linspace(0, 1.5e6, 151)[1:]*2*pi):
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
		
		if i1%10 == 0 and i2%10 == 0:
			plt.figure()
			plt.plot(t*1000, output.expect[1], label='Mechanical 1 num expec')
			plt.plot(t*1000, output.expect[2], label='Mechanical 2 num expec')
			plt.plot(t*1000, alpha_n1, label="Pulse 1")
			plt.plot(t*1000, alpha_n2, label="Pulse 2")
			plt.legend()
			plt.xlabel("Time (ms)")
			plt.title(r'$\omega_1/2\pi=\SI{%.3f}{MHz},\omega_2/2\pi=\SI{%.3f}{MHz}$' % (omega1*1e-6/(2*pi), omega2*1e-6/(2*pi)))
			plt.savefig(f'{plot_dir}/{i1}-{i2}.pdf', transparent=True, format='pdf', bbox_inches='tight')
			plt.close()

		omega1_vals.append(omega1)
		omega2_vals.append(omega2)
		n1_val = output.expect[1][-1]
		n2_val = output.expect[2][-1]
		n1_vals.append(n1_val)
		n2_vals.append(n2_val)
		with open(output_file, 'a') as f:
			f.write(f'{i1},{i2},{omega1},{omega2},{n1_val},{n2_val}\n')

x = array(omega1_vals).reshape(150,150) *1e-6/(2*pi)
y = array(omega2_vals).reshape(150,150) *1e-6/(2*pi)
z = array(n2_vals).reshape(150,150)

max_i = argmax(z)
max_x = x.flat[max_i]
max_y = y.flat[max_i]
print(r'Max <n2>: %.3f at %.3f MHz and %.3f MHz' % (z.flat[max_i], max_x, max_y))

fig, ax = plt.subplots()
mesh = ax.pcolormesh(x, y, z)
mesh.set_edgecolor('face') # https://stackoverflow.com/a/27096694
#plt.plot(max_x, max_y, 'ro')

ax.set_aspect('equal')
ax.set_xlabel(r'$\omega_1/2\pi$ (MHz)')
ax.set_ylabel(r'$\omega_2/2\pi$ (MHz)')

if os.path.exists('paper/omega1-omega2-sweep-zoomed.csv'):
	zoomed_data = loadtxt('paper/omega1-omega2-sweep-zoomed.csv', delimiter=',')
	zoomed_omega1_vals = list(zoomed_data[:, 2])
	zoomed_omega2_vals = list(zoomed_data[:, 3])
	zoomed_n1_vals = list(zoomed_data[:, 4])
	zoomed_n2_vals = list(zoomed_data[:, 5])
	x = array(zoomed_omega1_vals).reshape(200,200) *1e-6/(2*pi)
	y = array(zoomed_omega2_vals).reshape(200,200) *1e-6/(2*pi)
	z = array(zoomed_n2_vals).reshape(200,200)

	in_ax = ax.inset_axes([0.4, 0.1, 0.5, 0.5])
	in_ax.pcolormesh(x, y, z, cmap=mesh.cmap, norm=mesh.norm).set_edgecolor('face')
	for spine in in_ax.spines.values():
		spine.set_edgecolor("red")
	in_ax.set_aspect('equal')
	ax.indicate_inset_zoom(in_ax, edgecolor = 'red')

fig.colorbar(mesh, ax=ax).set_label(r'$\left<n_2\right>$')
plt.savefig(output_plot, transparent=True, format='pdf', bbox_inches='tight')
