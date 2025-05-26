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

data = np.loadtxt('output-approach-horizontal-1.8.csv', delimiter=',')
x = data[:,0] / (2*np.pi) / 1e6
y = data[:,1]

plt.plot(x, y)
plt.xlabel('$\omega_2/2\pi$ (MHz)')
plt.ylabel('Fidelity')
plt.show()
