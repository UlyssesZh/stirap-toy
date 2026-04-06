
import numpy as np
import scipy.constants as sc

def params_small():
    #N = dim(mechanical1) = dim(mechanical2), Nopt = dim(optical)
    N = 4
    Nopt = 2
    #single photon coupling strength of resonators
    g1 = 5 * np.pi
    g2 = 5 * np.pi
#20 #40(10K)
    #mechanical frequencies
    omega_1 = 1.2e6 * 2*np.pi
    omega_2 = 1.8e6 * 2*np.pi
    #was1.5e6
    #Cavity frequency, detector frequency
    omega_c = 5.4 * 2 * np.pi * 1e14

    #cavity and detector linewidths
#    kappa = np.pi * 0.5e5
    kappa = np.pi * 1.0e5
    kappa_d = np.pi * 1e4


    #Pump parameters: 40e-5 for stirap original
#    sigma_1 = 40e-5
#    sigma_2 =  40e-5
    #Pump parameters: 20e-5 for fstirap original
    sigma_1 = 60e-5
    sigma_2 =  60e-5
    sigma_blue = 20e-4
    sigma_readout = 5e-4
    #10k
#    sigma_1 = 40e-6
#    sigma_2 =  40e-6
    tau = sigma_1 / 1.43
    sigma_spring = 50e-3
#tau for stirap 1.43 (f1.3) (s1.6)
    #coherent amplitude of pumps (stirap)
#    alpha1 = 1000
#    alpha1 = 2000
#    alpha2 = 2000
    #coherent amplitude of pumps (fstirap)
    alpha1 = 2000
    alpha2 = 2000

    #pump frequencies, resonant (beamsplitter) case. Phi is the net phase in the exponenent of the operators
    n = 1#the coefficient in omega_c - n*omega_1
    omega_l1 = (omega_c - omega_1)
    omega_l2 = (omega_c - omega_2)

    omega_lspring = omega_c -n* omega_1
    omega_lspring2 = omega_c -n* omega_2

    phi_1 = omega_c - omega_l1 - omega_1
    phi_2 = omega_c - omega_l2 - omega_2
    phi_1_p = omega_c - omega_l1 + omega_1
    phi_2_p = omega_c - omega_l2 + omega_2
    phi_spring_1 = omega_c - omega_lspring - omega_1 #Note that this is the phi detuned from resonator 1
    phi_spring_2 = omega_c - omega_lspring - omega_2
    phi_spring_1_plus = (n+1) * omega_1
    phi_spring_1_minus = (n-1) * omega_1
    phi_spring_2_plus = (n*omega_1 +omega_2)
    phi_spring_2_minus = n*omega_1 - omega_2


    #Delta, according to the optomechanics bible
    delta = (omega_lspring - omega_c)
    delta2 = (omega_lspring2 - omega_c)

    #Environemnt temp, number of thermal phonons
    T = 1

    T1 = 40e-1
    T2 = 40e-3
    nth_1 = 1 / (np.exp((sc.hbar*omega_1) / (T * sc.k))-1)
    nth_2 = 1 / (np.exp((sc.hbar*omega_2) / (T * sc.k))-1)

    #thermal occupation at the mechanical frequency
    nb1 = sc.k * T1 / (sc.hbar * omega_1)
    nb2 = sc.k * T2 / (sc.hbar * omega_2)


    #Quality Factors
    Q1 = 1e9
    Q2 = 1e9

    #Membrane loss parameters
    gamma_1 = omega_1 / Q1
    gamma_2 = omega_2 / Q2


    sc.hbar = 1
    #phonon diffusion constant
    D1 = gamma_1 * nth_1
    D2 = gamma_2 * nth_2

    #effective masses, zero point motions
    meff1 = 1e-12
    meff2 = 1e-12
    x1zpf = np.sqrt(sc.hbar / (2 * meff1 * omega_1))
    x2zpf = np.sqrt(sc.hbar / (2 * meff2 * omega_2))

    #Optomechanical Swapping
    omega_bar = (omega_1 + omega_2) / 2
    delta_bar = 2.3e6 * 2 * np.pi
    J = 2 * g1 * g2 * np.sqrt(alpha1 * alpha2) * ((((omega_bar - delta_bar)/(((kappa**2))/4 + (omega_bar - delta_bar)**2)) - ((omega_bar + delta_bar)/((kappa**2))/4 + (omega_bar - delta_bar)**2)))


    tau_half_pi = 4*np.pi / (5*omega_1)
    #timing parameters orginal was 2e-3
#    t_delay = 0.008/2
    t_delay = 0
#time for STIRAP: oringal was 2e-3
#    t_max = 1e-3
#    t_min = -1e-3
#    t_min = -0.5e-3
#time used for readout
#    t_max = 5e-6
#    t_min = -2e-6
#    t_min = -1e-1
#    t_max = 1e-1
#time used for thermal equil
#    t_max = 1.11e1
#    t_max = 2e-2
#    t_min = -2e-2
    #timing parameters for mechanical swapping
#    t_max = 8e-6
#    t_min = -8e-6
    #num_step = 250000
#time used for blue detuning
#    t_min = -5e-3
#    t_max = 5e-3
#    t_min = -0.05e-3
#    t_max = 0.05e-3

#time for detection
#    t_max = 0.5e-7
#    t_min = -2e-7
#time for testing heating analytic results
#    t_max = 1.92e-4
#    t_min = 0.0
#time for fstirap
    t_max = 2.0e-3
    t_min = -2.0e-3
    num_step = 15001
    tlist  = np.linspace(t_min,t_max, num_step)

    params_dict = {
    "N": N,
    "Nopt": Nopt,
    "g1" : g1,
    "g2" : g2,
    "omega_1" : omega_1,
    "omega_2" : omega_2,
    "omega_c" : omega_c,
    "kappa" : kappa,
    "kappa_d" : kappa_d,
    "sigma_1" : sigma_1,
    "sigma_2" : sigma_2,
    "sigma_blue" : sigma_blue,
    "sigma_readout" : sigma_readout,
    "gamma_1" : gamma_1,
    "gamma_2" : gamma_2,
    "tau" : tau,
    "tau_half_pi" : tau_half_pi,
    "sigma_spring" : sigma_spring,
    "alpha1" : alpha1,
    "alpha2" : alpha2,
    "omega_l1" : omega_l1,
    "omega_l2" : omega_l2,
    "omega_lspring" : omega_lspring,
    "phi_1" : phi_1,
    "phi_2" : phi_2,
    "phi_1_p" : phi_1_p,
    "phi_2_p" : phi_2_p,
    "phi_spring_1" : phi_spring_1,
    "phi_spring_1+" : phi_spring_1_plus,
    "phi_spring_1-" : phi_spring_1_minus,
    "phi_spring_2+" : phi_spring_2_plus,
    "phi_spring_2-" : phi_spring_2_minus,
    "delta" : delta,
    "delta2" : delta2,
    "t_delay" : t_delay,
    "t_min" : t_min,
    "t_max" : t_max,
    "num_step" : num_step,
    "tlist" : tlist,
    "D1" : D1,
    "D2" : D2,
    "nth_1" : nth_1,
    "nth_2" : nth_2,
    "x1zpf" : x1zpf,
    "x2zpf" : x2zpf,
    "n" : n

    }
    return params_dict
