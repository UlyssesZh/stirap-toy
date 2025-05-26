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

data = np.load('main.npy')
x, y = np.meshgrid(np.linspace(1, 2, 201)+0.2, np.linspace(1, 2, 201)+0.8)
x = x.T
y = y.T
n1, n2 = np.meshgrid(np.arange(dim.n1), np.arange(dim.n2))

fig, axs = plt.subplots(dim.n1, dim.n2, figsize=(12, 9))
for n1 in range(dim.n1):
	for n2 in range(dim.n2):
		ax = axs[n1, n2]
		ax.set_title(r'$n_1 = %d, n_2 = %d$' % (n1, n2))
		mesh = ax.pcolormesh(x, y, data[:,:,n1,n2], vmin=0, vmax=1)
		if n2 == 0:
			ax.set_ylabel(r'$\omega_2/2\pi$ (MHz)')
		else:
			ax.get_yaxis().set_ticks([])
		if n1 == dim.n2-1:
			ax.set_xlabel(r'$\omega_1/2\pi$ (MHz)')
		else:
			ax.get_xaxis().set_ticks([])
		ax.set_aspect('equal', adjustable='box')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(mesh, cax=cbar_ax)

plt.show()
