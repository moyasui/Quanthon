import numpy as np
from itertools import product
from .mappers_utils import *
from .utils import get_all_paulis, is_hermitian, get_pauli



def jordan_wigner(hamiltonian):
    
    '''
    Perform the Jordan-Wigner transform for which maps second quantisation Hamiltonian to qubit Hamiltonians.
    N.B. 0th qubit to the left, always.
    arg:
        hamiltonian: an instance of the Hamiltonian class, the second quantisation Hamiltonian containing the overlap integrals.
    
    return:
        h_pauli: list of tuples, the qubit Hamiltonian in terms of Pauli strings.'''
    h_pauli = []
    n = hamiltonian.one_body_coeffs.shape[0]

    # one body
    for i in range(n):
        for j in range(n):

            if np.allclose(hamiltonian.one_body_coeffs[i, j],0):
                continue

            # diagonal terms
            if i == j:
                # II...III - I...IZ...I

                all_i = n * 'I'
                h_pauli.append((all_i, 0.5 * hamiltonian.one_body_coeffs[i, j]))

                z_term = list(all_i)
                # z_term[n-i-1] = 'Z'
                z_term[i] = 'Z'
                z_term = ''.join(z_term)

                h_pauli.append((z_term, -0.5 * hamiltonian.one_body_coeffs[i, j]))

            else:
                ob_op = get_ob_op(n, i, j, hamiltonian.one_body_coeffs[i, j])
                h_pauli.extend(ob_op)
     
                
    # two body
                
    # print("unhealthy terms after one body", _check_health_jw(h_pauli, n))
    # return h_pauli

    for p in range(n):
        for q in range(n):
            if p == q:
                continue
            for r in range(n):
                for s in range(n):
                    if r == s:
                        continue


                    coeff = hamiltonian.two_body_coeffs[p, q, r, s]

                    if np.allclose(coeff,0):
                        continue
                    
                    # print(f'pqrs: {p}{q}{r}{s}, coeff: {coeff}')

                    
                    # if the outer indices are not all greater than the inner indices
                    if (p > q and r > s) or (p < q and r < s):
                        coeff *= -1
                        # print(coeff)

                    if (p == r and q == s) or (p == s and q == r):

                        # print(p, q, r, s)
                        all_i = n*'I'

                        pz = list(n*'I')
                        pz[p] = 'Z'
                        pz = ''.join(pz)

                        qz = list(n*'I')
                        qz[q] = 'Z'
                        qz = ''.join(qz)

                        qpz = list(n*'I')
                        qpz[p] = 'Z'
                        qpz[q] = 'Z'
                        qpz = ''.join(qpz)

                        # print(p, q, r, s, coeff)
                        # print(all_i, pz, qz, qpz)
                        # print()
                        # print("all_i", all_i)
                        # print("pz", pz)
                        # print("qz", qz)
                        # print("qpz", qpz)

                        h_pauli.append((all_i, 0.25 * coeff))
                        h_pauli.append((pz, -0.25 * coeff))
                        h_pauli.append((qz, -0.25 * coeff))
                        h_pauli.append((qpz, 0.25 * coeff))

                        
                    elif len({p, q, r, s}) == 4:
                        # if they are all different

                        tb_op = get_tb_4u(n, p, q, r, s, coeff)
                        h_pauli.extend(tb_op)

                    elif len({p, q, r, s}) == 3:
                        # 0223
                        tb_op = get_tb_3u(n, p, q, r, s, coeff)
                        h_pauli.extend(tb_op)


    if check_health_jw(h_pauli, n):
        raise ValueError("There are terms whose length is not equal to the number of qubits.")
    
    h_pauli = simplify_pauli_terms(h_pauli)
    return h_pauli



def pauli_decomposition(h_mat):

    '''
    Perform the Pauli decomposition of a Hermitian matrix of size (2^n,2^n).
    args:
        h_mat: 2d array, the Hamiltonian matrix to be decomposed.
    
    return:
        h_pauli: list of tuples, the qubit Hamiltonian in terms of Pauli strings.
    
    ''' 

    # checks
    if not is_hermitian(h_mat):
        raise ValueError("The Hamiltonian matrix is not Hermitian.")
    
    h_pauli = []
    
    if not is_power_of_two(h_mat.shape[0]):
        h_mat = expand_h_mat(h_mat)
        
    n = int(np.log2(h_mat.shape[0]))
    all_paulis = get_all_paulis(n)
    
    for paulis in all_paulis:
        coeff = (1/2**n) * np.trace(get_pauli(paulis) @ h_mat)
        if coeff != 0:
            h_pauli.append((paulis, coeff))
    
    
    return h_pauli