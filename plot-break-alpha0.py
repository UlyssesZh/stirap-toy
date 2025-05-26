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

data = np.load('main-break-alpha0.npy')
x, y = np.meshgrid(np.linspace(1, 2, 51), np.linspace(1, 2, 51))
x = x.T
y = y.T
n1, n2 = np.meshgrid(np.arange(dim.n1), np.arange(dim.n2))
z1 = np.sum(data * n1, axis=(2,3))
z2 = np.sum(data * n2, axis=(2,3))

mesh = plt.pcolormesh(x, y, z1)
plt.colorbar(mesh).set_label(r'$\left<n_1\right>$')
plt.xlabel('$\omega_1/2\pi$ (MHz)')
plt.ylabel('$\omega_2/2\pi$ (MHz)')
plt.show()

mesh = plt.pcolormesh(x, y, z2)
plt.colorbar(mesh).set_label(r'$\left<n_2\right>$')
plt.xlabel('$\omega_1/2\pi$ (MHz)')
plt.ylabel('$\omega_2/2\pi$ (MHz)')
plt.show()
