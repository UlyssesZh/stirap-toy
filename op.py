from numpy import linspace, sin, cos, sqrt
from qutip import fock, fock_dm, tensor, create, destroy, qeye, qzero

from dim import *
from param import *
from t import *
from utils import expi, thermal, rev

a = tensor(destroy(nc), qeye(n1), qeye(n2))
ad = tensor(create(nc), qeye(n1), qeye(n2))
b1 = tensor(qeye(nc), destroy(n1), qeye(n2))
bd1 = tensor(qeye(nc), create(n1), qeye(n2))
b2 = tensor(qeye(nc), qeye(n1), destroy(n2))
bd2 = tensor(qeye(nc), qeye(n1), create(n2))
zero = tensor(qzero(nc), qeye(n1), qeye(n2))

H = [zero]
for oa, ca in zip((ad, a, ad, a), (alpha1*expi(Delta1*t), alpha1*expi(-Delta1*t), alpha2*expi(Delta2*t), alpha2*expi(-Delta2*t))):
	for ob, cb in zip((bd1, b1, bd2, b2), (g1*expi(omega1*t), g1*expi(-omega1*t), g2*expi(omega2*t), g2*expi(-omega2*t))):
		H.append([oa*ob, ca*cb])

c_ops=[
	sqrt(kappa)*a,
	sqrt(Gamma_m1*(1+nth1))*b1,
	sqrt(Gamma_m1*nth1)*bd1,
	sqrt(Gamma_m2*(1+nth2))*b2,
	sqrt(Gamma_m2*nth2)*bd2
]

e_ops=[
	ad*a,
	bd1*b1,
	bd2*b2
]

psi0 = tensor(fock(nc,0), fock(n1,1)-1j*fock(n1,0), fock(n2,0)).unit()
psi1 = tensor(fock_dm(nc,0), rev(psi0.ptrace(2)), rev(psi0.ptrace(1))).unit()
