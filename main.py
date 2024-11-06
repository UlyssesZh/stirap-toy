#!/usr/bin/env python

from numpy import linspace, pi, exp, sqrt, sin, cos
import matplotlib.pyplot as plt
from qutip import mesolve, fidelity, hinton, plot_wigner

from t import t, alpha_n1, alpha_n2
from op import H, c_ops, e_ops, psi0, psi1

import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

output = mesolve(H, psi0, t, c_ops, e_ops, options={
	#'nsteps': 2000000000,
	'progress_bar': 'tqdm',
	'store_states': True
})
final = output.states[-1]
print(f'Fidelity: {fidelity(psi1, final)}')

plt.plot(t, output.expect[1], label='Mechanical 1 num expec')
plt.plot(t, output.expect[2], label='Mechanical 2 num expec')
plt.plot(t, alpha_n1, label="Pulse 1")
plt.plot(t, alpha_n2, label="Pulse 2")
plt.legend()
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.show()

plot_wigner(final, linspace(-5,5,1000), linspace(-5,5,1000))
plt.show()
