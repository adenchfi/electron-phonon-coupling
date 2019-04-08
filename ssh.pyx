# this script is to test a simple brute force minimization of the free energy with respect to lattice bond lengths
# haven't changed anything from test_semiclassical_ssh.py yet
import networkx as nx
import numpy as np
cimport numpy as np
import sys
from scipy.sparse.linalg import eigs
import scipy.linalg as scilin
import numpy.linalg as lin
import matplotlib.pyplot as plt
import argparse
from latticeutil import create_inhomo_lattice
cimport cython
DTYPE = np.complex128
ctypedef np.complex128_t DTYPE_t

# evaluate action
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)   # Deactivate negative indexing.
cdef complex objective_eval(float beta, complex[:] myphiqs, complex[:] lat_disp, complex[:] elec_disp, complex[:,:,:] omega_1s, complex[:,:,:] omega_2s, int num, long[:, :] I, np.ndarray[DTYPE_t, ndim=2] xi): 
    # one doesn't have to explicitly construct the operator G, just the operator xi, since we've summed over Matsubara frequencies
    # compute Tr ln (G^{-1}_{el-ph})
    
    cdef complex elphterm = 0
    cdef complex j = 1j
    for eta in range(num):
        for etap in range(num):
            xi[eta, etap] = 0
            if (eta == etap):
                xi[eta, etap] = xi[eta, etap] + elec_disp[eta]
            elphterm = 0
            for gamma in range(num): ### TODO check if 1j multiplication belongs here - first indication is yes, it does
                elphterm = elphterm + j*((myphiqs[gamma].conjugate())*omega_1s[eta, etap, gamma] - myphiqs[gamma]*omega_2s[eta, etap, gamma])
            xi[eta, etap] = xi[eta, etap] - elphterm
            
    cdef complex trlnGinv = np.trace(scilin.logm(np.add(I, scilin.expm(-beta*xi))))  ### TODO - if no electrons, trlnGinv should reduce to ln(I) = 0, not ln(2I)
    cdef complex ph_term = 0
    for q in range(num):
        ph_term += (myphiqs[q].conjugate())*(lat_disp[q])*myphiqs[q]
    
    return (ph_term + trlnGinv) ### originally had minus, but maybe Matsubara sum changes it

@cython.boundscheck(False) # turn off bounds-checking for entire function
def start(args):

    dims = args.dims
    cdef int dimi = int(dims[0])
    cdef int dimj = int(dims[1])
    cdef int dimk = int(dims[2])
    disorder_type = args.disorder_type[0] # default disorder if I think it'll help convergence

    cdef float disorder_strength = float(args.disorder_strength[0])
    outname = args.outfile[0]
    cdef float alpha = float(args.alpha[0])
    cdef float t = float(args.t[0])
    cdef float K = float(args.K[0])
    cdef float C = float(args.C[0])
    cdef float Temp = float(args.Tem[0])
    choice = int(args.which[0])
    cdef float m = 1
    per = True
    if (int(args.periodic[0]) != 1):
        per = False
        

    # create graph
    inits, G = create_inhomo_lattice(dimi, dimj, dimk,None, 1, periodic=per)
    inits_w, G_w = create_inhomo_lattice(dimi, dimj, dimk,"Alternating", disorder_strength, periodic=per)

    #####
    # create graph Laplacian for lattice and find its eigenvectors and eigenvalues
    L_at_py = (K/2)*(nx.laplacian_matrix(G, weight='weight').real).toarray()
    cdef np.ndarray L_at = np.array(L_at_py, dtype=DTYPE)

    lat_eigvalspy, lat_eigvecspy = lin.eig(L_at)
    idxs = lat_eigvalspy.argsort()[::-1]
    lat_eigvalspy = np.array(lat_eigvalspy, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] lat_eigvals = lat_eigvalspy[idxs]
    cdef np.ndarray lat_eigvecs = (lat_eigvecspy[:, idxs])
    zero_check_idxs = np.abs(lat_eigvals) < 1e-14
    lat_eigvals[zero_check_idxs] = 0
    # create graph Laplacian for electrons; will need to add chemical potential as necessary, make other modifications
    H_TBpy = t*(np.real(nx.laplacian_matrix(G, weight=None))).toarray()
    cdef np.ndarray[DTYPE_t, ndim=2] H_TB = np.array(H_TBpy, dtype=DTYPE)
    el_energiespy, el_eigvecs = lin.eig(H_TB)
    idxs = el_energiespy.argsort()[::1]
    el_energiespy = np.array(el_energiespy, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] el_eigvals = el_energiespy[idxs]

    zero_check_idxs = np.abs(el_eigvals) < 1e-13
    el_eigvals[zero_check_idxs] = 0

    #####


    #####
    # initialize phiqs to be zero for now; may do otherwise when exploring metatsbale states
    #import numpy.random as rand
    cdef np.ndarray[DTYPE_t, ndim=1] phiqs = np.zeros((G.number_of_nodes()), dtype=DTYPE) ##
    cdef complex[:] phiqs_view = phiqs
    #phiqs = rand.randn(G.number_of_nodes(), 1)/50
    #phiqs = np.array(phiqs, dtype=complex)
    #####
    # initialize temperature and whatnot
    cdef float kb = 1

    cdef Py_ssize_t N = sum([1 for e in el_eigvals]) # one inefficient way of counting number of eigenvectors
    #####
    cdef Py_ssize_t M = (lat_eigvecs[:, 0]).size

    cdef complex[:] el_eigvals_view = el_eigvals
    cdef complex[:] lat_eigvals_view = lat_eigvals
    
    # compute omega function elements
    cdef np.ndarray[DTYPE_t, ndim=3] omega_1 = np.zeros((N, N, N), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=3] omega_2 = np.zeros((N, N, N), dtype=DTYPE)
    
    cdef complex[:,:,:] omega_1_view = omega_1
    cdef complex[:,:,:] omega_2_view = omega_2
    cdef complex[:,:] el_eigvecs_view = el_eigvecs
    cdef complex[:,:] lat_eigvecs_view = lat_eigvecs
    
    import time
    startt = time.clock()
    cdef complex eps_ijetaetap = 0
    # TODO - it appears these omegas are not generated correctly
    cdef Py_ssize_t eta, etap, i, j, gamma
    print("Beginning omega function generation")
    for eta in range(N):
        for etap in range(N):
            for i in range(M):
                for j in range(M):
                    eps_ijetaetap = ((el_eigvecs_view[i, eta]).conjugate())*el_eigvecs_view[j, etap]
                    for gamma in range(N):
                        omega_1_view[eta, etap, gamma] = omega_1_view[eta, etap, gamma] + alpha*((lat_eigvecs_view[i, gamma] - lat_eigvecs_view[j, gamma]).conjugate())*eps_ijetaetap

                        omega_2_view[eta, etap, gamma] = omega_2_view[eta, etap, gamma] + alpha*(lat_eigvecs_view[i, gamma] - lat_eigvecs_view[j, gamma])*eps_ijetaetap


    print("Time for omega function generation, ", time.clock() - startt)
    #####
    # Begin Loop
    cdef np.ndarray[DTYPE_t, ndim=2] xi = np.zeros((N,N), dtype=DTYPE)
    cdef np.ndarray I = np.eye(N, dtype=np.int)
    cdef long[:,:] I_view = I

    cdef complex curr_action = objective_eval(1/(kb*Temp), phiqs_view, lat_eigvals_view, el_eigvals_view, omega_1_view, omega_2_view, N, I_view, xi)
    cdef complex min_action = curr_action
    print("Starting Action:", curr_action)
    cdef np.ndarray min_phiqs = np.empty_like(phiqs, dtype=DTYPE)
    min_phiqs[:] = phiqs
    cdef complex[:] min_phiqs_view = min_phiqs
    print("Starting Phiqs:", phiqs)
    np.set_printoptions(linewidth=135)
    cdef double gam = 1
    cdef complex jay = 1j
    mingam = 1e-1
    
    for i in range(4): # 8 times we repeat this hackjob of a line search
        print("Pass: ", i)
        for idx1 in range(N): 
            gam = 0.2
            flag1 = False
            while (gam > mingam):
                if (flag1):
                    gam = gam/4
                    ####### we have a search direction; lat_eigvecs[:, idx1]. We would like to minimize our free energy along that direction before continuing in the next direction. 
                phiqs_view[idx1] = min_phiqs_view[idx1]

                phiqs_view[idx1] = phiqs_view[idx1] - gam*(jay)**(2*i+1) # we let it minimize across both real and imaginary parts of phi    
                curr_action = objective_eval(1/(kb*Temp), phiqs_view, lat_eigvals_view, el_eigvals_view, omega_1_view, omega_2_view, N, I_view, xi)

                if ((curr_action.real) - (min_action.real) < -1e-8):
                    print("Action", curr_action)
                    print("Accepted: ", idx1, " With Stepsize:", gam*(jay)**(2*i+1))

                    min_action = curr_action
                    min_phiqs_view[idx1] = phiqs_view[idx1]
                    flag1 = False
                else:
                    flag1 = True

    ##### End Loop

    min_phiqs[np.abs(min_phiqs) < 1e-14] = 0
    ##### Print Stuff
    print("Minimum action", min_action)
    print("Printing and saving Minimum phiqs setup")
    print(min_phiqs)
    name=("K%f-t%f-alpha%f-dims-%d-%d-%d-periodic-%d-Temp-%f" % (K, t, alpha, dimi, dimj, dimk, per, Temp))
    np.save(name, min_phiqs)
    #####


    ##### Compute average lattice displacements X_i
    # We can now compute the average lattice displacements from the above and compute the effective noninteracting electron DoS using a modified Hamiltonian

    ##### inverse graph Fourier transform to get real-space components 
    phirs = np.zeros((N,1), dtype=complex)
    for q, phi_q in enumerate(min_phiqs):
        np.add(phirs[:], np.reshape(min_phiqs[q]*lat_eigvecs[:, q], np.shape(phirs)), out=phirs[:])
        #####
        # Perhaps now that I have phi(l, tau=0) I must convert this to average lattice displacements. Maybe <x_i> = <phi_i e^{-phi_i}>?

    disps = np.zeros((N, 1), dtype=complex)
    for r, phi_r in enumerate(phirs):
        #disps[r] = phirs[r] * np.exp(np.conj(phirs[r])*phirs[r])  # weight
        ### TODO - not sure how to compute it, but need an imaginary unit??? supposed to be plus below, not minus
        disps[r] = 1j*(phirs[r] - (phirs[r]).conjugate())  # lattice displacement is sum of phi and phibar
    print("Real-space values of phi_l, disps")
    print(phirs)
    print(disps)

    ##### Modify TB Hamiltonian with lattice changes X_i - X_j
    H_TB = np.array(H_TB, dtype=complex)
    for i, rowi in enumerate(H_TB):
        for j, colj in enumerate(H_TB):
            if (i != j):
                if (H_TB[i, j] != 0):
                    if (i > j):
                        H_TB[i, j] = H_TB[i, j] + alpha*(disps[i] - disps[j])
                        H_TB[j, i]= H_TB[i, j]
                        #v = v[:, idx]

    ### testing purposes, do weighted laplacian
    w, v = lin.eig(H_TB)

    from quickdos import delta, dos
    oms = w
    #oms = np.sqrt(np.abs(oms))
    Es = np.linspace(min(oms) - 0.5, max(oms) + 0.5, 400)
    Es2 = np.linspace(min(w) - 0.5, max(w) + 0.5, 400)
    DOS = [dos(oms, E) for E in Es]
    DOS2 = [dos(w, E) for E in Es2]

    plt.figure(1)
    plt.plot(Es, DOS, label="Density of States of TB Ham")
    plt.title("Optimized SSH DOS", fontsize=24)
    plt.xlabel("Energy", fontsize=22)
    plt.ylim(bottom=0)
    plt.xticks(fontsize=18)
    plt.legend(fontsize=20)
    fign = name + ".png"
    plt.savefig(fign)

    #plt.figure(2)
    #plt.plot(Es2, DOS2, label="Density of States of Graph Laplacian")
    #plt.xlabel("Energy", fontsize=22)
    #plt.ylim(bottom=0)
    #plt.xticks(fontsize=18)
    #plt.legend(fontsize=20)


    plt.show()
    # want to get the eigenvectors of this graph and observe the (1) most common states (2) states with highest eigenvalue? 
