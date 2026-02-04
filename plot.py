#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

import dim
import param

plt.rcParams.update({
	'text.latex.preamble': r'\usepackage{lmodern}\usepackage{siunitx}',
	'text.usetex': True,
	'font.size': 11,
	'font.family': 'lmodern',
	'lines.linewidth': 1
})

'''
data = np.loadtxt('omega1-omega2-sweep_1', delimiter=',')
x = data[:,0].reshape(51,51) / (2*np.pi) / 1e6
y = data[:,1].reshape(51,51) / (2*np.pi) / 1e6
z = data[:,2].reshape(51,51)**2

max_i = np.argmax(z)
max_x = x.flat[max_i]
max_y = y.flat[max_i]
print('Max fidelity: %.2f at %.2f MHz and %.2f MHz' % (z.flat[max_i], max_x, max_y))

mesh = plt.pcolormesh(x, y, z)
plt.colorbar(mesh).set_label('Fidelity')
mesh.set_edgecolor('face') # https://stackoverflow.com/a/27096694
omega2 = param.omega2 / (2*np.pi) / 1e6
plt.plot([1.84, 1.92], [omega2, omega2], color='red')
plt.plot(max_x, max_y, 'ro')
plt.xlabel('$\omega_1/2\pi$ (MHz)')
plt.ylabel('$\omega_2/2\pi$ (MHz)')
plt.show()
'''

data = np.load('main.npy')
x, y = np.meshgrid(np.linspace(1, 2, 201)+0.2, np.linspace(1, 2, 201)+0.8)
x = x.T
y = y.T
n1, n2 = np.meshgrid(np.arange(dim.n1), np.arange(dim.n2))
z1 = np.sum(data * n1, axis=(2,3))
z2 = np.sum(data * n2, axis=(2,3))

mesh = plt.pcolormesh(x, y, z1)
plt.colorbar(mesh).set_label(r'$\left<n_1\right>$')
mesh.set_edgecolor('face') # https://stackoverflow.com/a/27096694
plt.xlabel('$\omega_1/2\pi$ (MHz)')
plt.ylabel('$\omega_2/2\pi$ (MHz)')
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

mesh = plt.pcolormesh(x, y, z2)
plt.colorbar(mesh).set_label(r'$\left<n_2\right>$')
mesh.set_edgecolor('face') # https://stackoverflow.com/a/27096694
plt.xlabel('$\omega_1/2\pi$ (MHz)')
plt.ylabel('$\omega_2/2\pi$ (MHz)')
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
