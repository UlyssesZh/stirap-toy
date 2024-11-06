from numpy import pi, exp, log, arange
from qutip import qdiags
from scipy.constants import hbar, k

def expi(x):
	return exp(1j*x)

def thermal_n(nbar, thres=0.01):
	return log(thres)/log(1-1/(1+nbar))

def thermal(nbar, n=thermal_n(nbar)):
	return qdiags((1-1/(1+nbar))**arange(n)).unit()

def rev(m):
	return 2*qdiags(m.diag()) - m

def bose(omega, T):
	return 1/(exp(hbar*omega/(k*T))-1)
