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

plot_dir = 'paper/sigma-omega-sweep'
output_file = 'paper/sigma-omega-sweep.csv'
output_plot = 'paper/sigma-omega-sweep.pdf'
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
	sigma_vals = []
	delta_omega_vals = []
	n1_vals = []
	n2_vals = []
else:
	sigma_vals = list(data[:, 2])
	delta_omega_vals = list(data[:, 3])
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
for i1, sigma in enumerate(linspace(0.2e-3, 1.2e-3, 51)):
	for i2, delta_omega in enumerate(linspace(0, 0.5e6, 51)*2*pi):
		count += 1
		if count <= len(data):
			continue

		param.sigma1 = param.sigma2 = sigma
		param.omega1 = 1.3e6*2*pi - delta_omega/2
		param.omega2 = 1.3e6*2*pi + delta_omega/2
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
			plt.title(r'$\sigma_1=\sigma_2=\SI{%.3f}{ms},(\omega_2-\omega_1)/2\pi=\SI{%.3f}{MHz}$' % (sigma*1e3, delta_omega*1e-6/(2*pi)))
			plt.savefig(f'{plot_dir}/{i1}-{i2}.pdf', transparent=True, format='pdf', bbox_inches='tight')
			plt.close()

		sigma_vals.append(sigma)
		delta_omega_vals.append(delta_omega)
		n1_val = output.expect[1][-1]
		n2_val = output.expect[2][-1]
		n1_vals.append(n1_val)
		n2_vals.append(n2_val)
		with open(output_file, 'a') as f:
			f.write(f'{i1},{i2},{sigma},{delta_omega},{n1_val},{n2_val}\n')

x = array(sigma_vals).reshape(51,51) *1e3
y = array(delta_omega_vals).reshape(51,51) *1e-6/(2*pi)
z = array(n2_vals).reshape(51,51)

max_i = argmax(z)
max_x = x.flat[max_i]
max_y = y.flat[max_i]
print(r'Max <n2>: %.3f at %.3f MHz and %.3f MHz' % (z.flat[max_i], max_x, max_y))

mesh = plt.pcolormesh(x, y, z)
plt.colorbar(mesh).set_label(r'$\left<n_2\right>$')
mesh.set_edgecolor('face') # https://stackoverflow.com/a/27096694
#plt.plot(max_x, max_y, 'ro')
plt.xlabel(r'$\sigma$ (ms)')
plt.ylabel(r'$(\omega_2-\omega_1)/2\pi$ (MHz)')
plt.savefig(output_plot, transparent=True, format='pdf', bbox_inches='tight')
