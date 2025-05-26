#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
	'text.latex.preamble': r'\usepackage{lmodern}\usepackage{siunitx}',
	'text.usetex': True,
	'font.size': 11,
	'font.family': 'lmodern',
	'lines.linewidth': 1
})

data = np.loadtxt('output-break-sigma.csv', delimiter=',')
x = data[:,0].reshape(51,51) / (2*np.pi) / 1e6
y = data[:,1].reshape(51,51) / (2*np.pi) / 1e6
z = data[:,2].reshape(51,51)**2

max_i = np.argmax(z)
max_x = x.flat[max_i]
max_y = y.flat[max_i]
print('Max fidelity: %.2f at %.2f MHz and %.2f MHz' % (z.flat[max_i], max_x, max_y))

mesh = plt.pcolormesh(x, y, z)
plt.colorbar(mesh).set_label('Fidelity')
plt.plot(max_x, max_y, 'ro')
plt.xlabel('$\omega_1/2\pi$ (MHz)')
plt.ylabel('$\omega_2/2\pi$ (MHz)')
plt.show()
