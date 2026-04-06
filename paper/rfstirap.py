

#qutip imports
from qutip import *
from qutip.ui.progressbar import BaseProgressBar, TextProgressBar
from qutip.tensor import flatten

#standard python imports
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.constants as sc
import sys as sys

#Our imports
from Decoherence import *
from pumps import *
#from params_small import *
from params_small_v import *
#from params_small_singlerot import *
from operators import *
from optical_spring_paramcheck import *
#from hamiltonian_optical_spring_pulse import *
#from hamiltonian_fstirap_reversedfstirap import *
#from hamiltonian_fstirap_reversedfstirap_triprot import *
#from hamiltonian_fstirap_reversedfstirap_constant import *
from hamiltonian_rfstirap import *
#from hamiltonian_stirap_nonrotating import *
#from hamiltonian_single_mode import *
np.set_printoptions(threshold=sys.maxsize)

def main():

    params_dict = params_small()
    #assign necessary paremeters from the dictionary
    N, Nopt, omega_1, omega_2, omega_c, omega_lspring, g1, g2, kappa, kappa_d, alpha1, alpha2, gamma_1, gamma_2,  sigma_1, sigma_2, tau, t_delay, t_min, t_max, numstep, tlist, D1, D2, nth_1, nth_2, x1zpf, x2zpf = params_dict["N"], params_dict["Nopt"], params_dict["omega_1"], params_dict["omega_2"],params_dict["omega_c"],params_dict["omega_lspring"], params_dict["g1"], params_dict["g2"], params_dict["kappa"],params_dict["kappa_d"], params_dict["alpha1"], params_dict["alpha2"],params_dict["gamma_1"],params_dict["gamma_2"], params_dict["sigma_1"], params_dict["sigma_2"], params_dict["tau"], params_dict["t_delay"], params_dict["t_min"], params_dict["t_max"], params_dict["num_step"], params_dict["tlist"], params_dict["D1"], params_dict["D2"], params_dict["nth_1"], params_dict["nth_2"], params_dict["x1zpf"], params_dict["x2zpf"]



    #Here we determine the analytic detuning parameters to be checked, note this is for pumping detuned from omega_1
    delta = omega_lspring - omega_c
    delta2 = (omega_c -omega_2 - omega_c)
    gamma_opt = get_gamma_opt(alpha1, g1, kappa, delta, omega_1)
    delta_omega = get_delta_omega(alpha1, g1, kappa, delta, omega_1)
    gamma1 = get_gamma_opt1(alpha1, g1, kappa, delta, omega_1)
    gamma2 = get_gamma_opt2(alpha2, g2, kappa, delta2, omega_2)
    delta1 = get_delta_omega1(alpha1, g1, kappa, delta, omega_1)
    delta2 = get_delta_omega2(alpha2, g2, kappa, delta, omega_2)
    print(nth_1 )
#    print(nth_1)
#    gamma_opt_trash = (0.82*alpha1*g1)**2 * ((kappa / ((kappa**2) / 4 + (delta + omega_1)**2)) - (kappa / ((kappa**2)/4+(delta - omega_1)**2)))
#    print(gamma_opt_trash)
    #gamma_opt_decay = plot_gamma_opt_decay(tlist, gamma_opt, tmin)
    #delta_omega_real = plot_real_spring_phase(tlist, delta_omega, tmin)

    #Get hamiltonian and operators from Hamiltonian fn, as well as figures
    fig, ax1 = plt.subplots(2,1)
    [H, a, b1, b2, d, fig, ax1] = Hamiltonian(params_dict, fig, ax1)
    p1_b1 = tensor(qeye(Nopt), basis(N, 1) * basis(N, 1).dag(), qeye(N))
    p1_b2 = tensor(qeye(Nopt), qeye(N), basis(N, 1) * basis(N, 1).dag())
#    print( a * (b1))
#    print(a.dag() * (b1.dag()))

    #testing coherent state plots
#    psic = coherent_dm(N,1)
#    hinton(psic)
#    plt.show()
#    print(a * (b1.dag()))


#    psi_0 =  0.5*(tensor(fock_dm(Nopt,0),fock_dm(N,1),fock_dm(N,0)) + tensor(fock_dm(Nopt,0),fock_dm(N,0),fock_dm(N,1)))
#    psi_0 =  (tensor(fock_dm(Nopt,0),fock_dm(N,1),fock_dm(N,0)) + tensor(fock_dm(Nopt,0),fock_dm(N,0),fock_dm(N,1)))
    psi_0 =  tensor(fock(Nopt,0),fock(N,1), fock(N,0)) + tensor(fock(Nopt,0),fock(N,0), fock(N,1))
    psi_0 =  ket2dm(psi_0)
#    print(psi_0)
#    psi_c = tensor(coherent_dm(Nopt,0),coherent_dm(N,1),coherent_dm(N,0))
#    psi_t = tensor(thermal_dm(Nopt,0),thermal_dm(N,3),thermal_dm(N,0))


#    psi_0 = tensor(fock(Nopt,0),fock(N,0),fock(N,0)) + tensor(fock(Nopt,0),fock(N,2),fock(N,0))
#    psi_0 = ket2dm(psi_0)
#    psi_0 = 0.5 * psi_0

#    psi_0 = tensor(fock_dm(Nopt,0), fock_dm(N,1), fock_dm(N,0))
#    psi_test = tensor(fock(Nopt,0),fock(N,0),fock(N,0)) - tensor(fock(Nopt,0),fock(N,0),fock(N,2))
#    psi_test = ket2dm(psi_test)
#    psi_test = 0.5 * psi_test
#    psi_test = tensor(fock_dm(Nopt,0), fock_dm(N,1), fock_dm(N,0))
    p = tensor(fock(N,0), fock(N,1)) - tensor(fock(N,1), fock(N,0))
    ps = tensor(fock(Nopt,0), p)
    psi_rf = 0.5 * ket2dm(ps)
    psi_plus = 0.5 * ((tensor(fock(Nopt,0),fock(N,0),fock(N,0)) + tensor(fock(Nopt,0),fock(N,1),fock(N,0))) + (tensor(fock(Nopt,0),fock(N,0),fock(N,0)) - tensor(fock(Nopt,0),fock(N,0),fock(N,1))))
#    psi_0 =  (1 / np.sqrt(2)) * ket2dm(((tensor(fock(Nopt,0),fock(N,0),fock(N,0)) + tensor(fock(Nopt,0),fock(N,1),fock(N,0))) + (tensor(fock(Nopt,0),fock(N,0),fock(N,0)) + tensor(fock(Nopt,0),fock(N,0),fock(N,1)))))
#    psi_0 = 0.5 * (tensor(fock(Nopt,0),fock(N,1),fock(N,0)) + (tensor(fock(Nopt,0),fock(N,0),fock(N,1))))

    psi_t = tensor(thermal_dm(Nopt, 0.001), thermal_dm(N, 0.5), thermal_dm(N, 0.5))
    psi_minus = qload('fstirap_result')

#    hinton(psi_0.ptrace(1))
#    plt.show()
#    hinton(psi_0.ptrace(2))
#    plt.show()

#    psi_test_c = tensor(coherent_dm(Nopt,0),coherent_dm(N,0),coherent_dm(N,3))
#    psi_test_t = tensor(thermal_dm(Nopt,0),thermal_dm(N,0),thermal_dm(N,3))

#    psi_1 = tensor(fock_dm(Nopt,0), fock_dm(N,1), fock_dm(N,0))
#    psi_none = tensor(fock_dm(Nopt,0), fock_dm(N,0), fock_dm(N,0))

    #red->blue->heralded density matrix
#    psi_h002 = 9.80586146e-01 * tensor(fock_dm(Nopt,0),fock_dm(N,1),fock_dm(N,0)) + 1.90373342e-2 * tensor(fock_dm(Nopt,0),fock_dm(N,2),fock_dm(N,0)) + 3.69217232e-04 * tensor(fock_dm(Nopt,0),fock_dm(N,3),fock_dm(N,0))
#    psi_h02 = 8.32042809e-01 * tensor(fock_dm(Nopt,0),fock_dm(N,1),fock_dm(N,0)) + 1.39783987e-01 * tensor(fock_dm(Nopt,0),fock_dm(N,2),fock_dm(N,0)) + 2.34570810e-02 * tensor(fock_dm(Nopt,0),fock_dm(N,3),fock_dm(N,0)) + 3.93187567e-03 * tensor(fock_dm(Nopt,0),fock_dm(N,4),fock_dm(N,0)) + 9.02596129e-01 * tensor(fock_dm(Nopt,0),fock_dm(N,0),fock_dm(N,0)) + 8.79650413e-02 * tensor(fock_dm(Nopt,0),fock_dm(N,0),fock_dm(N,1))
#    psi_h02_true = psi_h02.unit()
#    psi_h02 = tensor(fock_dm(Nopt,0), 8.31398514e-01 * fock_dm(N,1)+ 1.40283588e-01 * fock_dm(N,2) + 2.35861204e-02 * fock_dm(N,3) + 3.95190235e-03 * fock_dm(N,4), 8.93106023e-01 * fock_dm(N,0) + 9.54730029e-02 * fock_dm(N,1)+1.02063447e-02 * fock_dm(N,2))
##    psi_h02 = psi_h02 = tensor(fock_dm(Nopt,0), (2/17)*(8.31398514e-01 * fock_dm(N,0)+ 1.40283588e-01 * fock_dm(N,1) + 2.35861204e-02*fock_dm(N,2))+(15/17)*(8.31398514e-01 * fock_dm(N,1)+ 1.40283588e-01 * fock_dm(N,2) + 2.35861204e-02 * fock_dm(N,3) + 3.95190235e-03 * fock_dm(N,4)), 8.93106023e-01 * fock_dm(N,0) + 9.54730029e-02 * fock_dm(N,1)+1.02063447e-02 * fock_dm(N,2))
#    psi_h02_true = psi_h02.unit()
#    print(psi_h02_true)
#    psi_h02 = 8.32042809e-01 * tensor(fock(Nopt,0),fock(N,1),fock(N,0)) + 1.39783987e-01 * tensor(fock(Nopt,0),fock(N,2),fock(N,0)) + 2.34570810e-02 * tensor(fock(Nopt,0),fock(N,3),fock(N,0)) + 3.93187567e-03 * tensor(fock(Nopt,0),fock(N,4),fock(N,0)) + 9.02596129e-01 * tensor(fock(Nopt,0),fock(N,0),fock(N,1)) + 8.79650413e-02 * tensor(fock(Nopt,0),fock(N,0),fock(N,2))
#    psi_h02 =  ket2dm(psi_h02)
#    psi_h1 = 4.98794484e-01 * tensor(fock_dm(Nopt,0),fock_dm(N,1),fock_dm(N,0)) + 2.50644153e-01 * tensor(fock_dm(Nopt,0),fock_dm(N,2),fock_dm(N,0)) + 1.25623381e-01 * tensor(fock_dm(Nopt,0),fock_dm(N,3),fock_dm(N,0)) + 6.28000589e-02 * tensor(fock_dm(Nopt,0),fock_dm(N,4),fock_dm(N,0)) + 3.13130983e-02 * tensor(fock_dm(Nopt,0),fock_dm(N,5),fock_dm(N,0)) + 1.55728542e-02 * tensor(fock_dm(Nopt,0),fock_dm(N,6),fock_dm(N,0)) + 7.72478467e-03 * tensor(fock_dm(Nopt,0),fock_dm(N,7),fock_dm(N,0)) + 3.82190928e-03 * tensor(fock_dm(Nopt,0),fock_dm(N,8),fock_dm(N,0))
#    print(psi_0)
#    print(psi_1)
#    print(psi_test)
    #increase the step number of the solver
    options = Options()
    options.nsteps = 2000000000
    nskips = 1000
#    print(nth_1)
#    print(nth_2)
#    print(D2)

#    rho_m1_initial=psi_h02.ptrace(1)
#    print(rho_m1_initial)
#    rho_m2_initial=psi_h02.ptrace(2)
#    print(rho_m2_initial)
    #Loop through detunings of the stokes pump

    #find which times to plot, on a frequency scale relevant to the delta_m
    #delta_t_prime = 2 * np.pi * 0.1 / np.absdelta_omega
    #nskips = int(round(delta_t_prime / ((t_max - t_min)/num_step)))
    #print(nskips)


        #omega_l2 = omega_l2 + (omega_c - omega_2)*0.1
    #Set to reuse hamiltonian data, and solve

    #Decay/thermalisation Operators
    cc = (np.sqrt(kappa)/5)* a
    cm1 = np.sqrt(gamma_1) * b1
    cm2 = np.sqrt(gamma_2) * b2
    c_thermal1 = (b1 + b1.dag()) * np.sqrt(D1)
    c_thermal2 = (b2 + b2.dag()) * np.sqrt(D2)
    c_true_thermal1a = np.sqrt(D1) * b1
    c_true_thermal1b = np.sqrt(gamma_1) * b1
    c_true_thermal1c = np.sqrt(D1) * b1.dag()
    c_true_thermal2a = np.sqrt(D2) * b2
    c_true_thermal2b = np.sqrt(gamma_2) * b2
    c_true_thermal2c = np.sqrt(D2) * b2.dag()


#    c_op_list =[]
#    c_op_list =[cc]
    c_op_list = [cc, c_thermal1, c_thermal2]
#    c_op_list = [cc, c_true_thermal1a, c_true_thermal1b, c_true_thermal1c, c_true_thermal2a, c_true_thermal2b, c_true_thermal2c]
#    c_op_list = [c_true_thermal1a, c_true_thermal1c, c_true_thermal2a, c_true_thermal2c]
#    c_op_list = [c_true_thermal1a, c_true_thermal1b, c_true_thermal1c, c_true_thermal2a, c_true_thermal2b, c_true_thermal2c]

    opts = Options(rhs_reuse=True)
    result = mesolve(H, psi_minus , tlist, c_op_list, [], options=options, progress_bar = True)
    qsave(result.states[-1], 'rfstirap_result')
    hinton(result.states[0].ptrace([1,2]))
#    result_nint = mesolve(a.dag()*a-a.dag()*a, psi_h02, tlist, c_op_list, [], options=options, progress_bar = True)
    fid = fidelity(tensor(fock_dm(Nopt,0), fock_dm(N,1), fock_dm(N,0)), result.states[-1])

    coherence_l1 = [np.sum(np.abs(rho.data.toarray())) - np.trace(rho.data.toarray()).real
                for rho in result.states]
    idx_10 = N  # Cavity=0, Mech1=1, Mech2=0
    idx_01 = 1  # Cavity=0, Mech1=0, Mech2=1

    # Extract the off-diagonal element modulus |rho_ab| at each time step
    # Note: state.full() converts the Qobj to a numpy array for fast indexing
    coherence_modulus = [np.abs(state[idx_10, idx_01]) for state in result.states]
    np.savez('coherence_evolution_rfstirap.npz',
         time=tlist,
         coherence=coherence_l1,
         params={'g1': g1, 'kappa': kappa}) # Optional: save metadata
    # Plotting the evolution
    plt.figure(figsize=(8, 5))
    plt.plot(tlist, coherence_l1, label=r'$|\rho_{(010),(001)}|$ (Coherence)', color='tab:orange', lw=2)

    # Optional: Plot populations for context
    pop_10 = [np.abs(state[idx_10, idx_10]) for state in result.states]
    pop_01 = [np.abs(state[idx_01, idx_01]) for state in result.states]
    plt.plot(tlist, pop_10, '--', label='Population |0,1,0>', alpha=0.6)

    plt.xlabel('Time')
    plt.ylabel('Modulus Value')
    plt.title('Evolution of Mechanical Coherence')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

        # --- 1. Compute Negativity ---
    # We ptrace indices 1 and 2 (the mechanical modes)
    # Then calculate negativity between index 0 and 1 of that reduced system
    def calculate_manual_negativity(rho):
        # 1. Isolate the two mechanical modes
        rho_mech = rho.ptrace([1, 2])

        # 2. Transpose the first resonator (index 0)
        # mask [1, 0] means transpose the first, keep the second
        rho_pt = partial_transpose(rho_mech, [1, 0])

        # 3. Get eigenvalues
        evals = rho_pt.eigenenergies()

        # 4. Sum only the absolute values of the negative eigenvalues
        # This is the formal definition of Negativity
        return np.sum(np.abs(evals[evals < 0]))

    # Apply to your 15,001 steps
    neg_evolution = [calculate_manual_negativity(rho) for rho in result.states]

    # --- 2. Extract Populations for context ---
    # Using expectation operators is faster, but since you have the states in memory:
    # We'll use the kets to ensure we are picking the right Fock states
    s1 = ket([0, 1, 0], [Nopt, N, N])
    s2 = ket([0, 0, 1], [Nopt, N, N])

    pop_10 = expect(s1 * s1.dag(), result.states)
    pop_01 = expect(s2 * s2.dag(), result.states)

    # --- 3. Save Data ---
    np.savez('negativity_evolution_rfstirap.npz',
             time=tlist,
             negativity=neg_evolution,
             pop_10=pop_10,
             pop_01=pop_01)

    # --- 4. Plotting ---
    plt.figure(figsize=(8, 5))

    # Plot Negativity (Black Dashed Line)
    plt.plot(tlist, neg_evolution, label='Negativity $\mathcal{N}$', color='black', linestyle='--', lw=2)

    # Plot Populations for context
    plt.plot(tlist, pop_10, label='Population $|0,1,0\\rangle$', color='tab:red', alpha=0.7)
    plt.plot(tlist, pop_01, label='Population $|0,0,1\\rangle$', color='tab:blue', alpha=0.7)

    plt.xlabel('Time (ms)')
    plt.ylabel('Value')
    plt.title('Evolution of Entanglement (Negativity) at 1K')
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


    rho_m1_final=result.states[-1].ptrace(1)
    print(rho_m1_final)
    rho_m2_final=result.states[-1].ptrace(2)
    print(rho_m2_final)
    rho_c_final=result.states[-1].ptrace(0)
    print(rho_c_final)
#    print(fidelity(fock_dm(N,1), 0.5*result.states[-1].ptrace(1)))

#    print("Expectation Value of Resonator 1=", expect(b1.dag() * b1,result.states[0]))
#    print("Expectation Value of Resonator 2=", expect(b2.dag() * b2,result.states[0]))

    state100 = []
    state200 = []
    state111 = []
    state211 = []
    for i in range(0,15001):
#        state10 = result.states[i].ptrace(1)[0,0]
#        state100.append(state10)
#        state20 = result.states[i].ptrace(2)[0,0]
#        state200.append(state20)
        state11 = result.states[i].ptrace(1)[1,1]
        state111.append(state11)
        state21 = result.states[i].ptrace(2)[1,1]
        state211.append(state21)







    print(fid)
    ax1[0].legend(loc = 'upper right')
    #ax3.legend(loc = 'upper center')
    #ax1[1].legend(loc = 'upper right')

    #ax3.tick_params(axis='y', colors = 'orange')
    #ax2.set_ylabel(r'$g\alpha (Hz)$', fontsize=20)
    #ax2.legend(loc = 'upper center')

    #lbls_list = [['$\ket{{}}$'.append(n) for n in range(N)], ["u", "d"]]

    #ax1[0,1].plot(hinton(result.states[-1]))
    plotlist=np.linspace(-2.0,2.0, 15001)
    Data_csv = np.column_stack((plotlist,expect(p1_b1.dag() * p1_b1, result.states),expect(p1_b2.dag() * p1_b2, result.states),fSTIRAP_P(-tlist, tau, t_delay, sigma_1),fSTIRAP_S(-tlist, tau, t_delay, sigma_2)))
    file_data_store('rfstirap_0+1_10mk.csv', Data_csv, sep=',')
    #print(expect(b2.dag() * b2, result.states[-1]))
    #plt.table(cellText = [])
    def plot_pop(result):
        #use 2 when generating actual stirap
        fig, ax = plt.subplots()
        ax.plot(plotlist, expect(b1.dag() * b1, result.states), marker = '', color = 'red',  linestyle = 'solid', label = 'Mechanical Resonator 1')
        ax.plot(plotlist, expect(b2.dag() * b2, result.states), marker = '', color = 'blue',  linestyle = 'solid', label = 'Mechanical Resonator 2')
#        ax.plot(plotlist, state111, marker = '', color = 'red',  linestyle = 'solid', label = 'Mode 1')
#        ax.plot(plotlist, state211, marker = '', color = 'blue',  linestyle = 'solid', label = 'Mode 2')
#        ax.plot(plotlist, expect(a.dag() * a, result.states), marker = '', color = 'green',  linestyle = 'dashed', label = 'Optical Cavity')
        ax.plot(plotlist, fSTIRAP_P(-tlist, tau, t_delay, sigma_1), marker = '', color = 'red',  linestyle = 'dotted',label = 'P Pulse')
        ax.plot(plotlist, fSTIRAP_S(-tlist, tau, t_delay, sigma_2), marker = '', color = 'blue',  linestyle = 'dotted', label = 'S Pulse')
#        ax.plot(tlist, state100, marker = '', color = 'pink',  linestyle = 'dashed', label = 'Mode 1')
#        ax.plot(plotlist, constant_pulse(tlist)-state100, marker = '', color = 'red',  linestyle = 'dotted', label = 'Mode 1')
#        ax.plot(plotlist, constant_pulse(tlist)-state200, marker = '', color = 'blue',  linestyle = 'dotted', label = 'Mode 2')
#        ax.plot(plotlist, state111, marker = '', color = 'red',  linestyle = 'solid', label = 'Mode 1')
#        ax.plot(plotlist, state211, marker = '', color = 'blue',  linestyle = 'solid', label = 'Mode 2')
#        ax.plot(plotlist, gaussian1(plotlist, 1000*tau, 1000*t_delay, 1000*sigma_1), marker = '', color = 'red',label = 'S Pulse')
#        ax.plot(plotlist, gaussian2(plotlist, 1000*tau, 1000*t_delay, 1000*sigma_1), marker = '', color = 'blue', label = 'P Pulse')
#        ax.plot(tlist, state1fid, marker = '', color = 'green',  linestyle = 'solid', label = 'fidelity')
#        ax.plot(plotlist, state1fid_nint, marker = '', color = 'red',  linestyle = 'dashed', label = 'fidelity')
#        ax.plot(plotlist, state2fid_nint, marker = '', color = 'blue',  linestyle = 'dashed', label = 'fidelity')
#        ax.plot(tlist, 1*constant_pulse(tlist)-statec00, marker = '', color = 'green',  linestyle = 'dashed', label = 'Cavity')
#        ax.plot(tlist, statec11, marker = '', color = 'green',  linestyle = 'solid', label = 'Cavity')
#        ax.errorbar(tlist, constant_pulse(tlist)-state1, yerr = asymmetric_error, xerr = None, ecolor= 'lightgrey', lolims=0, marker = '', color = 'Red',  linestyle = 'dashed', label = 'Mode 1')
#        ax.errorbar(tlist, constant_pulse(tlist)-state2, yerr = asymmetric_error, xerr = None, ecolor= 'lightblue',lolims=0, marker = '', color = 'Blue',  linestyle = 'dashed', label = 'Mode 2')
#        ax.plot(tlist, gaussian1(tlist, tau, t_delay, sigma_1), marker = '', color = 'red',label = 'Constant Pulse')
#        ax.plot(tlist,gaussian2(tlist, tau, t_delay, sigma_1), marker = '', color = 'blue', label = 'P Pulse')
#        ax.plot(tlist, (1 - np.exp(-1*D2*tlist)), marker = '', color = 'blue', label = '$e^{-n_{th} \gamma_{m}t}$')
#        ax.plot(tlist, (fSTIRAP_P(-tlist,  tau, t_delay, sigma_1)+fSTIRAP_P(tlist,  tau, t_delay, sigma_1)), marker = '', color = 'red',label = 'P Pulse')
#        ax.plot(tlist, (fSTIRAP_S(-tlist,  tau, t_delay, sigma_2)+fSTIRAP_S(tlist,  tau, t_delay, sigma_2)), marker = '', color = 'blue', label = 'S Pulse')
#        ax.legend(loc = 'upper right')
#        ax.set_ylim(bottom=0.0, top=1.01)
#        ax.set_xlim(left=-0.5, right=0.5)
        ax.set_xlabel('Time(ms)', fontsize = 14)
        ax.set_ylabel('Probability Amplitude', fontsize = 14)
#        ax.set_title('Reverse fStirap |01><01| + |10><10|', fontsize = 20)
        ax.set_title('', fontsize = 20)
        plt.savefig('C:/Users/Ian Hedgepeth/Documents/STIRAP plots/Plots for supplemental/100821/rfSTIRAP_100821_kappa_nRWA_(-).eps', format='eps')
        plt.show()
    plot_pop(result)

    hinton(result.states[0])
    plt.show()
    hinton(result.states[-1])
    plt.show()
    hinton(result.states[-1].ptrace(1))
    plt.show()
    hinton(result.states[-1].ptrace(2))
    plt.show()






if __name__ == '__main__':
    main()
