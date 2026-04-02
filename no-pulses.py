#!/usr/bin/env python

from numpy import linspace, loadtxt, array
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

plot_dir = 'paper'

param.kappa = 0
H, c_ops, e_ops, psi0, psi1 = setup()
e_ops = [
	qutip.tensor(qutip.qeye(dim.nc), qutip.basis(dim.n1,1) * qutip.basis(dim.n1,1).dag(), qutip.qeye(dim.n2)),
	qutip.tensor(qutip.qeye(dim.nc), qutip.qeye(dim.n1), qutip.basis(dim.n2,1) * qutip.basis(dim.n2,1).dag()),
]

options={
	#'nsteps': 2000000000,
	'progress_bar': 'tqdm',
	'store_final_state': True
}

t, alpha_n1, alpha_n2 = t_dep.t, t_dep.alpha_n1, t_dep.alpha_n2
output = mesolve(H, psi0, t, c_ops, e_ops, options=options)
final = output.final_state
fid = fidelity(psi1, final)
print(f'n1_val = {output.expect[0][-1]}')
print(f'n2_val = {output.expect[1][-1]}')
print(f'Fidelity: {fid**2}')

plt.figure()
plt.plot(t*1000, output.expect[0], label='Mechanical 1 prob')
plt.plot(t*1000, output.expect[1], label='Mechanical 2 prob')
plt.plot(t*1000, alpha_n1, label="Pulse 1")
plt.plot(t*1000, alpha_n2, label="Pulse 2")
plt.ylim(0, 2)
plt.xlabel("Time (ms)")
plt.ylabel("Amplitude")
plt.savefig(f'{plot_dir}/no-pulses.pdf', transparent=True, format='pdf', bbox_inches='tight')
plt.close()
