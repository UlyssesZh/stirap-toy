

#qutip imports
from qutip import *
from qutip.ui.progressbar import BaseProgressBar, TextProgressBar
from qutip.tensor import flatten


#standard python imports
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import numpy as np
import scipy.constants as sc
import sys as sys
import pandas as pd
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.ticker import FixedFormatter, FixedLocator

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
from hamiltonian_fstirap import *
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


    # Define the basis states for a 2-mode system
    # fock(N, n) creates a state with n particles in a Hilbert space of size N
    state_10 = tensor(fock(Nopt,0), fock(N, 1), fock(N, 0)) # |1,0>
    state_01 = tensor(fock(Nopt,0), fock(N, 0), fock(N, 1)) # |0,1>

    # Create the Bell state (10 + 01) / sqrt(2)
    psi_bell = (state_10 - state_01).unit()

    # Calculate the density matrix rho = |psi><psi|
    rho_bell = ket2dm(psi_bell)

#    psi_0 =  (tensor(fock_dm(Nopt,0),fock_dm(N,1),fock_dm(N,0)) + tensor(fock_dm(Nopt,0),fock_dm(N,0),fock_dm(N,1)))
#    psi_0 =  0.5*(tensor(fock_dm(Nopt,0),fock_dm(N,1),fock_dm(N,0)) + tensor(fock_dm(Nopt,0),fock_dm(N,0),fock_dm(N,0)))
    psi_0 =  1/np.sqrt(2) *(tensor(fock(Nopt,0),fock(N,0), fock(N,0)) + tensor(fock(Nopt,0),fock(N,1), fock(N,0)))

    psi_0 =  ket2dm(psi_0)
    psi_0 =  tensor(fock_dm(Nopt,0),fock_dm(N,1),fock_dm(N,0))
#    print(psi_0)
    psi_c = tensor(coherent_dm(Nopt,0),coherent_dm(N,1),coherent_dm(N,0))
    psi_c_2 = tensor(coherent_dm(Nopt,0),coherent_dm(N,0.707),coherent_dm(N,-0.707))

    psi_t = tensor(thermal_dm(Nopt,0),thermal_dm(N,1),thermal_dm(N,0))

#    psi_0 = tensor(fock(Nopt,0),fock(N,0),fock(N,0)) + tensor(fock(Nopt,0),fock(N,1),fock(N,0))
#    psi_0 = ket2dm(psi_0)
#    psi_0 =  psi_0
#    psi_0 = tensor(fock_dm(Nopt,0),fock_dm(N,0),fock_dm(N,0)) + tensor(fock_dm(Nopt,0),fock_dm(N,1),fock_dm(N,0))
#    psi_0 = tensor(maximally_mixed_dm(Nopt),maximally_mixed_dm(N),maximally_mixed_dm(N))
##    psi_0 = np.sqrt(0.5)*(tensor(fock(Nopt,0),fock(N,0),fock(N,0)) - tensor(fock(Nopt,0),fock(N,1),fock(N,0)))
#    psi_0 = tensor(fock(Nopt,0),fock(N,5),fock(N,0))

#    psi_0 = 0.5 * psi_0

#    psi_0 = tensor(fock_dm(Nopt,0), fock_dm(N,1), fock_dm(N,0))
#    psi_test = tensor(fock(Nopt,0),fock(N,0),fock(N,0)) - tensor(fock(Nopt,0),fock(N,0),fock(N,2))
#    psi_test = ket2dm(psi_test)
#    psi_test = 0.5 * psi_test
#    psi_test = tensor(fock_dm(Nopt,0), fock_dm(N,1), fock_dm(N,0))
#    psi_test_c = tensor(coherent_dm(Nopt,0),coherent_dm(N,0),coherent_dm(N,3))
#    psi_test_t = tensor(thermal_dm(Nopt,0),thermal_dm(N,0),thermal_dm(N,3))

#    psi_1 = tensor(fock_dm(Nopt,0), fock_dm(N,1), fock_dm(N,0))
#    psi_none = tensor(fock_dm(Nopt,0), fock_dm(N,0), fock_dm(N,0))

    #red->blue->heralded density matrix
#    psi_h002 = 9.80586146e-01 * tensor(fock_dm(Nopt,0),fock_dm(N,1),fock_dm(N,0)) + 1.90373342e-2 * tensor(fock_dm(Nopt,0),fock_dm(N,2),fock_dm(N,0)) + 3.69217232e-04 * tensor(fock_dm(Nopt,0),fock_dm(N,3),fock_dm(N,0))
#    psi_h02_unnorm = 8.32042809e-01 * tensor(fock_dm(Nopt,0),fock_dm(N,1),fock_dm(N,0)) + 1.39783987e-01 * tensor(fock_dm(Nopt,0),fock_dm(N,2),fock_dm(N,0)) + 2.34570810e-02 * tensor(fock_dm(Nopt,0),fock_dm(N,3),fock_dm(N,0)) + 3.93187567e-03 * tensor(fock_dm(Nopt,0),fock_dm(N,4),fock_dm(N,0))
#    psi_h02 = psi_h02_unnorm / psi_h02_unnorm.tr()
#    psi_h02_true = psi_h02.unit()
#    psi_h02 = tensor(fock_dm(Nopt,0), 8.31398514e-01 * fock_dm(N,1)+ 1.40283588e-01 * fock_dm(N,2) + 2.35861204e-02 * fock_dm(N,3) + 3.95190235e-03 * fock_dm(N,4), 8.93106023e-01 * fock_dm(N,0) + 9.54730029e-02 * fock_dm(N,1)+1.02063447e-02 * fock_dm(N,2))
##    psi_h02 = psi_h02 = tensor(fock_dm(Nopt,0), (2/17)*(8.31398514e-01 * fock_dm(N,0)+ 1.40283588e-01 * fock_dm(N,1) + 2.35861204e-02*fock_dm(N,2))+(15/17)*(8.31398514e-01 * fock_dm(N,1)+ 1.40283588e-01 * fock_dm(N,2) + 2.35861204e-02 * fock_dm(N,3) + 3.95190235e-03 * fock_dm(N,4)), 8.93106023e-01 * fock_dm(N,0) + 9.54730029e-02 * fock_dm(N,1)+1.02063447e-02 * fock_dm(N,2))
#    psi_h02_true = psi_h02.unit()
#    print(psi_h02_true)
#    psi_h02 = 8.32042809e-01 * tensor(fock(Nopt,0),fock(N,1),fock(N,0)) + 1.39783987e-01 * tensor(fock(Nopt,0),fock(N,2),fock(N,0)) + 2.34570810e-02 * tensor(fock(Nopt,0),fock(N,3),fock(N,0)) + 3.93187567e-03 * tensor(fock(Nopt,0),fock(N,4),fock(N,0)) + 9.02596129e-01 * tensor(fock(Nopt,0),fock(N,0),fock(N,1)) + 8.79650413e-02 * tensor(fock(Nopt,0),fock(N,0),fock(N,2))
#    psi_h02 =  ket2dm(psi_h02)
#    psi_h1 = 4.98794484e-01 * tensor(fock_dm(Nopt,0),fock_dm(N,1),fock_dm(N,0)) + 2.50644153e-01 * tensor(fock_dm(Nopt,0),fock_dm(N,2),fock_dm(N,0)) + 1.25623381e-01 * tensor(fock_dm(Nopt,0),fock_dm(N,3),fock_dm(N,0)) + 6.28000589e-02 * tensor(fock_dm(Nopt,0),fock_dm(N,4),fock_dm(N,0)) + 3.13130983e-02 * tensor(fock_dm(Nopt,0),fock_dm(N,5),fock_dm(N,0)) + 1.55728542e-02 * tensor(fock_dm(Nopt,0),fock_dm(N,6),fock_dm(N,0)) + 7.72478467e-03 * tensor(fock_dm(Nopt,0),fock_dm(N,7),fock_dm(N,0)) + 3.82190928e-03 * tensor(fock_dm(Nopt,0),fock_dm(N,8),fock_dm(N,0))


    p = tensor(fock(N,0), fock(N,1)) + tensor(fock(N,1), fock(N,0))
    ps = tensor(fock(Nopt,0), p)
    p = tensor(fock(N,0), fock(N,1), fock(N,0))+  tensor(fock(N,0), fock(N,0), fock(N,0))
    psi_rf = 0.5 * ket2dm(p)
#    hinton(psi_rf)

    #increase the step number of the solver
    options = Options()
    options.nsteps = 2000000000
    nskips = 1000


    #Decay/thermalisation Operators
    cc = (np.sqrt(kappa)/5)* a
    cm1 = np.sqrt(gamma_1) * b1
    cm2 = np.sqrt(gamma_2) * b2
    c_thermal1_1 = b1 * np.sqrt(D1+gamma_1)
    c_thermal1_2 = b1.dag() * np.sqrt(D1)
    c_thermal2_1 = b2 * np.sqrt(D2+gamma_2)
    c_thermal2_2 = b2.dag() * np.sqrt(D2)
    c_min_cooling_1 = np.sqrt(1) * b1.dag()
    c_min_cooling_2 = np.sqrt(1) * b2.dag()
    c_true_thermal1a = np.sqrt(D1) * b1
    c_true_thermal1b = np.sqrt(gamma_1) * b1
    c_true_thermal1c = np.sqrt(D1) * b1.dag()
    c_true_thermal2a = np.sqrt(D2) * b2
    c_true_thermal2b = np.sqrt(gamma_2) * b2
    c_true_thermal2c = np.sqrt(D2) * b2.dag()


    c_op_list =[]
#    c_op_list =[cc]
    c_op_list = [cc, c_thermal1_1, c_thermal1_2, c_thermal2_1, c_thermal2_2]
#    c_op_list = [cc, c_true_thermal1a, c_true_thermal1b, c_true_thermal1c, c_true_thermal2a, c_true_thermal2b, c_true_thermal2c]
#    c_op_list = [c_true_thermal1a, c_true_thermal1c, c_true_thermal2a, c_true_thermal2c]
#    c_op_list = [c_true_thermal1a, c_true_thermal1b, c_true_thermal1c, c_true_thermal2a, c_true_thermal2b, c_true_thermal2c]

    opts = Options(rhs_reuse=True)
    result = mesolve(H, psi_0 , tlist, c_op_list, [], options=options, progress_bar = True)

##
    # Manually check one state at the peak of your transfer
#    rho_mid = result.states[len(result.states)//2]
#    rho_mech = rho_mid.ptrace([1, 2])

    # Perform partial transpose on the first mechanical mode (index 0)
#    rho_pt = partial_transpose(rho_mech, [1, 0])

    # Look at the eigenvalues
#    evals = rho_pt.eigenenergies()
#    print("Eigenvalues of Partial Transpose:", evals)
##


    qsave(result.states[-1], 'fstirap_result')
    hinton(result.states[0].ptrace([1,2]))
    idx_10 = N  # Cavity=0, Mech1=1, Mech2=0
    idx_01 = 1  # Cavity=0, Mech1=0, Mech2=1
    # Assuming 'states' is the list of density matrices from mesolve
    # We use a lambda to sum the absolute values of off-diagonals: sum|rho_ij| for i != j
    coherence_l1 = [np.sum(np.abs(rho.data.toarray())) - np.trace(rho.data.toarray()).real
                for rho in result.states]
    # Extract the off-diagonal element modulus |rho_ab| at each time step
    # Note: state.full() converts the Qobj to a numpy array for fast indexing
    coherence_modulus = [np.abs(state[idx_10, idx_01]) for state in result.states]
    np.savez('coherence_evolution_fstirap.npz',
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
    np.savez('negativity_evolution_fstirap.npz',
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
    #Calculate the Fidelity
    fid = fidelity(rho_bell, result.states[-1])
    print(fid)


#    rho_m1_final=result.states[-1].ptrace(1)
#    print(rho_m1_final)
#    rho_m2_final=result.states[-1].ptrace(2)
#    print(rho_m2_final)
#    rho_c_final=result.states[-1].ptrace(0)
#    rho_m1_m2_final= ptrace(result.states[-1],[1,2])
#    print(rho_m1_m2_final)
#    rho_final=result.states[-1]
#    print(rho_final)

    file_data_store('final_state_qutip.csv',result.states[-1].full() , sep=',')



#    state100 = []
#    state200 = []
#    state111 = []
#    state211 = []
#    for i in range(0,15000):
#        state10 = result.states[i].ptrace(1)[0,0]
#        state100.append(state10)
#        state20 = result.states[i].ptrace(2)[0,0]
#        state200.append(state20)
#        state11 = result.states[i].ptrace(1)[1,1]
    #    state111.append(state11)
#        state21 = result.states[i].ptrace(2)[1,1]
#        state211.append(state21)






    # 1. Define local operators for an N-dimensional mode
    # We use the projection into the {0,1} subspace
    def local_op(theta):
        # Projections
        p0 = fock(N, 0) * fock(N, 0).dag()
        p1 = fock(N, 1) * fock(N, 1).dag()

        # Pauli Z equivalent: |0><0| - |1><1|
        sz = p0 - p1
        # Pauli X equivalent: |0><1| + |1><0|
        # This can also be written as (a + a.dag()) projected to this subspace
        sx = fock(N, 0) * fock(N, 1).dag() + fock(N, 1) * fock(N, 0).dag()

        return np.cos(theta) * sz + np.sin(theta) * sx

    # 2. Embed local measurements into the 3-mode Hilbert space [Nopt, N1, N2]
    # Alice measures N1 (mode index 1), Bob measures N2 (mode index 2)
    A1 = tensor(qeye(Nopt), local_op(0), qeye(N))
    A2 = tensor(qeye(Nopt), local_op(np.pi/2), qeye(N))

    B1 = tensor(qeye(Nopt), qeye(N), local_op(np.pi/4))
    B2 = tensor(qeye(Nopt), qeye(N), local_op(3*np.pi/4))

    # 3. Construct the CHSH Bell Operator
    # S = E(a,b) - E(a,b') + E(a',b) + E(a',b')
    CHSH_op = A1*B1 - A1*B2 + A2*B1 + A2*B2

    # 4. Calculate S
    S_value = expect(CHSH_op, result.states[-1])
    print(f"CHSH Value S: {S_value.real:.4f}")




    xvec = np.linspace(-2.5, 2.5, 150)
    Wigner = wigner(result.states[-1], xvec, xvec)
    wmap = wigner_cmap(Wigner)
    nrm = mpl.colors.Normalize(-Wigner.max(), Wigner.max())

    ax1[0].legend(loc = 'upper right')


    plotlist  = np.linspace(-2.0,2.0, 15001)
    Data_csv = np.column_stack((plotlist,expect(p1_b1.dag() * p1_b1, result.states),expect(p1_b2.dag() * p1_b2, result.states),fSTIRAP_P(-tlist, tau, t_delay, sigma_1),fSTIRAP_S(-tlist, tau, t_delay, sigma_2)))
    file_data_store('fstirap_0+1_10mk.csv', Data_csv, sep=',')

    def plot_pop(result):
        #use 2 when generating actual stirap
        fig= plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)
        ax2 = ax.twinx()
#        ax.plot(plotlist, expect(b1.dag() * b1, result.states), marker = '', color = 'red',  linestyle = 'solid', label = 'Mechanical Mode 1')
#        ax.plot(plotlist, expect(b2.dag() * b2, result.states), marker = '', color = 'blue',  linestyle = 'solid', label = 'Mechanical Mode 2')
        ax.plot(plotlist, expect(p1_b1, result.states), marker = '', color = 'red',  linestyle = 'solid', label = 'Mechanical Mode 1')
        ax.plot(plotlist, expect(p1_b2, result.states), marker = '', color = 'blue',  linestyle = 'solid', label = 'Mechanical Mode 2')

#        ax.plot(plotlist, expect(a.dag() * a, result.states), marker = '', color = 'green',  linestyle = 'solid', label = 'Optical Cavity')
        ax2.plot(plotlist, 2000*fSTIRAP_P(tlist, tau, t_delay, sigma_1), marker = '', color = 'red',  linestyle = 'dotted',label = 'P Pulse')
        ax2.plot(plotlist, 2000*fSTIRAP_S(tlist, tau, t_delay, sigma_2), marker = '', color = 'blue',  linestyle = 'dotted', label = 'S Pulse')
#        ax.plot(tlist, state100, marker = '', color = 'pink',  linestyle = 'dashed', label = 'Mode 1')
#        ax.plot(plotlist, constant_pulse(tlist)-state100, marker = '', color = 'red',  linestyle = 'dotted', label = 'Mode 1')
#        ax.plot(plotlist, constant_pulse(tlist)-state200, marker = '', color = 'blue',  linestyle = 'dotted', label = 'Mode 2')
#        ax.plot(plotlist, state111, marker = '', color = 'red',  linestyle = 'solid', label = 'Mode 1')
#        ax.plot(plotlist, state211, marker = '', color = 'blue',  linestyle = 'solid', label = 'Mode 2')
#        ax.plot(plotlist, gaussian1(plotlist, 1000*tau, 1000*t_delay, 1000*sigma_1), marker = '', color = 'red', linestyle = 'dotted', label = 'S Pulse')
#        ax.plot(plotlist, gaussian2(plotlist, 1000*tau, 1000*t_delay, 1000*sigma_1), marker = '', color = 'blue', linestyle = 'dotted', label = 'P Pulse')
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
        ax.set_xlabel('Time(ms)', fontsize = 20)
        ax.set_ylabel('Probability Amplitude', fontsize = 20)
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax2.set_ylabel('Pulse Amplitude', fontsize=20)
        ax2.tick_params(axis='y', labelsize=18)
        ax2.set_ylim(bottom=0.0, top=4000)
        ax.set_title('', fontsize = 25)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.savefig('C:/Users/Ian Hedgepeth/Documents/atom/stirap/STIRAP/V3 plots/fstirap_1K_0+1.pdf', format='pdf', bbox_inches='tight')
        plt.show()
    plot_pop(result)
    def format_hinton(fig, ax, filename):
        # 1. Manually fix the color scaling for the squares (PolyCollection)
        mappable = None
        for art in ax.get_children():
            # Hinton plots use PolyCollection; we set the color limits here
            if hasattr(art, 'set_clim'):
                art.set_clim(-0.6, 0.6)
                mappable = art

        # 2. Set spatial axis limits
        ax.set_xlim([3.8, 7.2])
        ax.set_ylim([18, 21.2])

        # 3. Set font sizes for the state labels [cite: 18, 20]
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(40)

        # 4. FIXED: Colorbar Alignment and Labels
        if len(fig.axes) > 1:
            cb_ax = fig.axes[-1]

            # Define the exact tick positions and correct the signs
            ticks = [-1.0, 0, 1.0]
            labels = ["-1.0", "0", "1.0"] # Manually fixing the top tick sign

            cb_ax.yaxis.set_major_locator(FixedLocator(ticks))
            cb_ax.yaxis.set_major_formatter(FixedFormatter(labels))

            # Ensure the visual axis of the colorbar matches the data range
            cb_ax.set_ylim(-0.6, 0.6)
            cb_ax.tick_params(labelsize=40)

        fig.savefig(filename, format='pdf', bbox_inches='tight')
    # Process the states (remove the vmin/vmax from the hinton call)
    fig1, ax1 = hinton(result.states[-1].ptrace([1,2]))
    format_hinton(fig1, ax1, 'C:/Users/Ian Hedgepeth/Documents/atom/stirap/STIRAP/V3 plots/initial_fSTIRAP_10mK_density.pdf')

    plt.show()

    def plot_Wigner(result):
        fig, ax = plt.subplots(figsize=(10, 8)) # Adjusted size for larger fonts

        # Calculate Wigner function
        W = wigner(result.states[-1], xvec, xvec)

        # 1. Establish symmetric limits for the colorbar
        limit = max(abs(W.min()), abs(W.max()))

        # 2. Use a Diverging Norm to force 0 to the center of the colormap
        norm = colors.TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)

        # Using 'RdBu_r' for a standard diverging colormap (Red positive, Blue negative)
        cmap = 'RdBu_r'

        # 3. Use imshow with interpolation for a smooth variation
        im = ax.imshow(
            W,
            cmap=cmap,
            norm=norm,
            interpolation='bilinear',
            origin='lower',
            extent=[xvec.min(), xvec.max(), xvec.min(), xvec.max()]
        )

        # 4. Formatting Ticks
        ax.tick_params(axis='both', which='major', labelsize=40)

        # 5. Symmetrical Colorbar with precise zero-anchoring
        cb1 = fig.colorbar(im, ax=ax, ticks=[-limit, 0, limit])

        # Format labels for readability at large sizes
        cb1.ax.set_yticklabels([f'{-limit:.2f}', '0', f'{limit:.2f}'])
        cb1.ax.tick_params(labelsize=40)

        # Remove colorbar "seams" for PDF export
        cb1.solids.set_edgecolor("face")

        fig.tight_layout()

        # Save and show
        save_path = 'C:/Users/Ian Hedgepeth/Documents/atom/stirap/STIRAP/V3 plots/Wigner_Plot_single_fstirap.pdf'
        plt.savefig(save_path, format='pdf', bbox_inches='tight')
        plt.show()


    plot_Wigner(result)


    Wigner_bell = wigner(bell_state(state='11'), xvec, xvec)
    wmap_bell = wigner_cmap(Wigner_bell)
    nrm_bell = mpl.colors.Normalize(-Wigner_bell.max(), Wigner_bell.max())
    psi_bell=bell_state(state='11')

#    def plot_Wigner_bell(psi_bell):
#        fig, ax = plt.subplots()
#        plt1 = ax.contourf(xvec, xvec, Wigner_bell, 100, cmap=wmap_bell)
#        ax.set_title("");
#        cb1 = fig.colorbar(plt1)
#        fig.tight_layout()
#        plt.show()
#    plot_Wigner_bell(psi_bell)

    def plot_Wigner_difference(result, psi_bell):
        # Calculate Wigner functions
        # Assuming xvec is defined in the global scope based on your snippet
        W1 = wigner(result.states[-1], xvec, xvec)
        W2 = wigner(psi_bell, xvec, xvec)
        Wigner_diff = W1 - W2

        fig, ax = plt.subplots(figsize=(10, 8))

        # 1. Establish symmetric limits for the colorbar
        limit = max(abs(W1.min()), abs(W1.max()))

        # 2. Use a Diverging Norm to force 0 to the center of the colormap
        norm = colors.TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)

        # Using 'RdBu_r' for a standard diverging colormap (Red positive, Blue negative)
        cmap = 'RdBu_r'

        # 3. Use imshow with interpolation for a smooth variation
        im = ax.imshow(
            Wigner_diff,
            cmap=cmap,
            norm=norm,
            interpolation='bilinear',
            origin='lower',
            extent=[xvec.min(), xvec.max(), xvec.min(), xvec.max()]
        )

        # 4. Formatting Ticks
        ax.tick_params(axis='both', which='major', labelsize=40)

        # 5. Symmetrical Colorbar with precise zero-anchoring
        cb1 = fig.colorbar(im, ax=ax, ticks=[-limit, 0, limit])

        # Format labels for readability at large sizes
        cb1.ax.set_yticklabels([f'{-limit:.2f}', '0', f'{limit:.2f}'])
        cb1.ax.tick_params(labelsize=40)

        # Remove colorbar "seams" for PDF export
        cb1.solids.set_edgecolor("face")

        fig.tight_layout()

        # Save and show
        save_path = 'C:/Users/Ian Hedgepeth/Documents/atom/stirap/STIRAP/V3 plots/Wigner_Difference_Smooth.pdf'
        plt.savefig(save_path, format='pdf', bbox_inches='tight')
        plt.show()
    plot_Wigner_difference(result, psi_bell)

if __name__ == '__main__':
    main()
