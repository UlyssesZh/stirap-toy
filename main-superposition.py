#!/usr/bin/env python

from numpy import linspace, loadtxt, array
import matplotlib.pyplot as plt
from qutip import mesolve, fidelity, hinton, plot_wigner, fock, tensor

import dim
import param
import t_dep
from op import setup

import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

#plt.rcParams.update({
#	'text.latex.preamble': r'\usepackage{lmodern}\usepackage{siunitx}',
#	'text.usetex': True,
#	'font.size': 11,
#	'font.family': 'lmodern',
#	'lines.linewidth': 1
#})


H, c_ops, e_ops, initial1, expected1 = setup()
initial1 = tensor(fock(dim.nc,0), fock(dim.n1,0)+fock(dim.n1,1), fock(dim.n2,0)).unit()
expected1 = tensor(fock(dim.nc,0), fock(dim.n1,0), fock(dim.n2,0)-fock(dim.n2,1)).unit()
initial2 = tensor(fock(dim.nc,0), fock(dim.n1,0)-fock(dim.n1,1), fock(dim.n2,0)).unit()
expected2 = tensor(fock(dim.nc,0), fock(dim.n1,0), fock(dim.n2,0)+fock(dim.n2,1)).unit()

options={
	#'nsteps': 2000000000,
	'progress_bar': 'tqdm',
	'store_final_state': True
}

t, alpha_n1, alpha_n2 = t_dep.t, t_dep.alpha_n1, t_dep.alpha_n2
output = mesolve(H, initial1, t, c_ops, e_ops=e_ops, options=options)
final1 = output.final_state
fid = fidelity(expected1, final1)
print(f'Fidelity: {fid**2}')

plt.plot(t*1000, output.expect[1], label='Mechanical 1 num expec')
plt.plot(t*1000, output.expect[2], label='Mechanical 2 num expec')
plt.plot(t*1000, alpha_n1, label="Pulse 1")
plt.plot(t*1000, alpha_n2, label="Pulse 2")
plt.ylim(0, 2)
plt.xlabel("Time (ms)")
plt.ylabel("Amplitude")
plt.title(r'$(|0\rangle+|1\rangle)/\sqrt{2}$')
plt.show()

output = mesolve(H, initial2, t, c_ops, e_ops=e_ops, options=options)
final2 = output.final_state
fid = fidelity(expected2, final2)
print(f'Fidelity: {fid**2}')

plt.plot(t*1000, output.expect[1], label='Mechanical 1 num expec')
plt.plot(t*1000, output.expect[2], label='Mechanical 2 num expec')
plt.plot(t*1000, alpha_n1, label="Pulse 1")
plt.plot(t*1000, alpha_n2, label="Pulse 2")
plt.ylim(0, 2)
plt.xlabel("Time (ms)")
plt.ylabel("Amplitude")
plt.title(r'$(|0\rangle-|1\rangle)/\sqrt{2}$')
plt.show()
