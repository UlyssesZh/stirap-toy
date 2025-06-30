#!/usr/bin/env python

from numpy import linspace, pi, exp, sqrt, sin, cos, around, ix_
import matplotlib.pyplot as plt
from qutip import mesolve, fidelity, hinton, plot_wigner, tensor, fock, fock_dm, ket2dm
from qutip.core.dimensions import to_tensor_rep

import dim
import param
import t_dep
from op import setup
from utils import rev

import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

plt.rcParams.update({
	'text.latex.preamble': r'\usepackage{lmodern}',
	'text.usetex': True,
	'font.size': 11,
	'font.family': 'lmodern',
	'lines.linewidth': 1
})

nc, n1, n2 = dim.nc, dim.n1, dim.n2
k00 = tensor(fock(nc,0), fock(n1,0), fock(n2,0))
k01 = tensor(fock(nc,0), fock(n1,0), fock(n2,1))
k10 = tensor(fock(nc,0), fock(n1,1), fock(n2,0))
psi0_list = [
	k10.unit(),
	(k00+k10).unit(),
	(k00-k10).unit(),
	(k00+1j*k10).unit(),
]
psi1_list = [
	(k10-k01).unit(),
	(sqrt(2)*k00+k10-k01).unit(),
	(-sqrt(2)*k00+k10-k01).unit(),
	(-1j*sqrt(2)*k00+k10-k01).unit(),
]
title_list = [
	r'$|10\rangle$',
	r'$|00\rangle+|10\rangle$',
	r'$|00\rangle-|10\rangle$',
	r'$|00\rangle+i|10\rangle$',
]
param.theta = pi/4

actual_psi1_list = []
for i, (psi0, psi1, title) in enumerate(zip(psi0_list, psi1_list, title_list)):
	H, c_ops, e_ops, *_ = setup()

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
	actual_psi1_list.append(final)

	plt.figure()
	plt.plot(t*1000, output.expect[1], label='Mechanical 1 num expec')
	plt.plot(t*1000, output.expect[2], label='Mechanical 2 num expec')
	plt.plot(t*1000, alpha_n1, label="Pulse 1")
	plt.plot(t*1000, alpha_n2, label="Pulse 2")
	plt.legend()
	plt.xlabel("Time (ms)")
	plt.ylabel("Amplitude")
	plt.title(title)
	plt.savefig(f'plot/{i}.pdf', transparent=True, format='pdf', bbox_inches='tight')
	plt.close()

param.tau1 = -param.tau1
for i, (psi0, psi1, title) in enumerate(zip(psi0_list, actual_psi1_list, title_list)):
	H, c_ops, e_ops, *_ = setup()

	options={
		#'nsteps': 2000000000,
		'progress_bar': 'tqdm',
		'store_final_state': True
	}

	t, alpha_n1, alpha_n2 = t_dep.t, t_dep.alpha_n1, t_dep.alpha_n2
	output = mesolve(H, psi1, t, c_ops, e_ops=e_ops, options=options)
	final = output.final_state
	fid = fidelity(psi0, final)
	print(f'Fidelity: {fid**2}')

	plt.figure()
	plt.plot(t*1000, output.expect[1], label='Mechanical 1 num expec')
	plt.plot(t*1000, output.expect[2], label='Mechanical 2 num expec')
	plt.plot(t*1000, alpha_n1, label="Pulse 1")
	plt.plot(t*1000, alpha_n2, label="Pulse 2")
	plt.legend()
	plt.xlabel("Time (ms)")
	plt.ylabel("Amplitude")
	plt.title(title)
	plt.savefig(f'plot/{i}-reversed.pdf', transparent=True, format='pdf', bbox_inches='tight')
	plt.close()
