#!/usr/bin/env python

from numpy import linspace, pi, exp, sqrt, sin, cos, loadtxt, array, argmax, arcsinh, log
import matplotlib.pyplot as plt
from qutip import mesolve, fidelity, hinton, plot_wigner
from scipy.special import lambertw

import dim
import param
import t_dep
from op import setup

import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

plot_dir = 'paper/tau-sigma-sweep-fractional'
output_file = 'paper/tau-sigma-sweep-fractional.csv'
output_plot = 'paper/tau-sigma-sweep-fractional.pdf'
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
	tau_vals = []
	n1_vals = []
	n2_vals = []
else:
	sigma_vals = list(data[:, 2])
	tau_vals = list(data[:, 3])
	n1_vals = list(data[:, 4])
	n2_vals = list(data[:, 5])

plt.rcParams.update({
	'text.latex.preamble': r'\usepackage{lmodern}\usepackage{siunitx}',
	'text.usetex': True,
	'font.size': 11,
	'font.family': 'lmodern',
	'lines.linewidth': 1
})

param.theta = pi/4
count = 0
for i1, sigma in enumerate(linspace(2e-4, 12e-4, 51)):
	for i2, tau in enumerate(linspace(2e-4, 12e-4, 51)):
		count += 1
		if count <= len(data):
			continue

		param.sigma1 = param.sigma2 = sigma
		param.tau1 = tau
		param.tau2 = -tau
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
			plt.title(r'$\sigma_1=\sigma_2=\SI{%.3f}{ms},\tau=\SI{%.3f}{ms}$' % (sigma*1000, tau*1000))
			plt.savefig(f'{plot_dir}/{i1}-{i2}.pdf', transparent=True, format='pdf', bbox_inches='tight')
			plt.close()

		sigma_vals.append(sigma)
		tau_vals.append(tau)
		n1_val = output.expect[1][-1]
		n2_val = output.expect[2][-1]
		n1_vals.append(n1_val)
		n2_vals.append(n2_val)
		with open(output_file, 'a') as f:
			f.write(f'{i1},{i2},{sigma},{tau},{n1_val},{n2_val}\n')

x = array(sigma_vals).reshape(51,51) *1000
y = array(tau_vals).reshape(51,51) *1000
z = array(n2_vals).reshape(51,51)

max_i = argmax(z)
max_x = x.flat[max_i]
max_y = y.flat[max_i]
print(r'Max <n2>: %.2f at %.2f ms and %.2f ms' % (z.flat[max_i], max_x, max_y))

fig, ax = plt.subplots()

mesh = ax.pcolormesh(x, y, z)
plt.colorbar(mesh).set_label(r'$\left<n_2\right>$')
mesh.set_edgecolor('face') # https://stackoverflow.com/a/27096694
#ax.plot(max_x, max_y, 'ro')

x1, x2 = ax.get_xlim()
y1, y2 = ax.get_ylim()

Omega0 = param.alpha0 * param.g1 * 2
n = 5
sigma = linspace(1e-4, 13e-4, 121)
theta = param.theta
tau_lo = (sqrt(log(2) + 2*arcsinh(cos(theta/2))) - sqrt(log(2))) * sigma/2
tau_hi = sqrt(2*lambertw(cos(theta/2)**2 / (n*sin(theta/2)) * Omega0 * sigma)) * sigma/2
ax.plot(sigma*1e3, tau_lo*1e3, 'blue')
ax.plot(sigma*1e3, tau_hi*1e3, 'blue')

ax.set_xlim(x1, x2)
ax.set_ylim(y1, y2)

ax.set_xlabel(r'$\sigma$ (ms)')
ax.set_ylabel(r'$\tau$ (ms)')
fig.savefig(output_plot, transparent=True, format='pdf', bbox_inches='tight')
