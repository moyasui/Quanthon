import numpy as np


def get_ob_op(n, i, j, coeff):
    # when the indices are different

    all_op = []

    if i < j:
        base = ['' if k == i or k == j else 'Z' if i < k < j else 'I' for k in range(n)]
        pauli_terms = [('XX', 0.25), ('XY', 0.25j), ('YX', -0.25j), ('YY', 0.25)]
    else:
        i, j = j, i  # the qubits that's smaller always on the left!!
        base = ['' if k == i or k == j else 'Z' if i < k < j else 'I' for k in range(n)] 
        # can't move outside conditional because i, j are swapped.
        pauli_terms = [('XX', 0.25), ('XY', -0.25j), ('YX', 0.25j), ('YY', 0.25)]

        # all gates acting on smaller numbers always on the left

    for op, factor in pauli_terms:
        new_op = base.copy()
        new_op[i] = op[0]
        new_op[j] = op[1]
        new_cf = factor * coeff
        all_op.append((''.join(new_op), new_cf))
    return all_op


def gen_prod_4u(ops):
    """
    Generate products of creation and annihilation operators in the given order.

    args:
        ops: list of length 4, consists of creation and annihilation operators depending on the order 
        of their indices.
            
    return:
        list of strings, each string is a term of the product
    """

    op1, op2, op3, op4 = ops
    products = []

    for t1, factor_1 in op1:
        for t2, factor_2 in op2:
            for t3, factor_3 in op3:
                for t4, factor_4 in op4:
                    op = t1 + t2 + t3 + t4 
                    factor = factor_1 * factor_2 * factor_3 * factor_4
                    products.append((op, factor))

    return products

def gen_prod_3u(ops):
    """
    Generate products of creation and annihilation operators in the given order.

    args:
        ops: list of length 3, consists of creation and annihilation operators depending on the order 
        of their indices.
            
    return:
        list of strings, each string is a term of the product
    """

    op1, op2, op3 = ops
    products = []

    for t1, factor_1 in op1:
        for t2, factor_2 in op2:
            for t3, factor_3 in op3:
                op = t1 + t2 + t3
                factor = factor_1 * factor_2 * factor_3
                products.append((op, factor))

    return products

def get_tb_base(n,ordered_indices):

    # get the I's and the Z's for the qubits not acted on. 
    tb_op_base = list(n * 'I')

    for i in range(n):
        if i in ordered_indices:
            continue

        if ordered_indices[0] < i < ordered_indices[1]:
            tb_op_base[i] = 'Z'
        
        if ordered_indices[2] < i < ordered_indices[3]:
            tb_op_base[i] = 'Z'
        
    return tb_op_base
    
def get_tb_4u(n,p,q,r,s,coeff):

    '''Get the two body operaors when all 4 indices are unique.'''
    
    ordered_indices = sorted([p,q,r,s])

    tb_op_base = get_tb_base(n, ordered_indices)

    creation = [('X', 0.5), ('Y', -0.5j)]
    annihilation = [('X', 0.5), ('Y', 0.5j)]

    sorted_zip = sorted(zip([p,q,r,s], [creation, creation, annihilation, annihilation]))
    ops = [op for _, op in sorted_zip]

    op_products = gen_prod_4u(ops)
    # print(op_products)

    all_tb_op = []
    for op, factor in op_products:
        # every product turn into one operator in the sum
        if np.allclose(coeff, 0):
            continue # skip coeffs close to 0
        
        new_op = tb_op_base.copy()
        
        for i, index in enumerate(ordered_indices):
            new_op[index] = op[i]
        
        new_op = ''.join(new_op)
        all_tb_op.append((new_op, factor * coeff))
        # print(new_op, coeff)
    return all_tb_op


def get_tb_3u(n, p, q, r, s, coeff):

    '''Get the two body operaors when 3 indices are unique.
    args:
        n: int, number of qubits,
        p, q, r, s: int, the indices of the qubits the second quantisation operation is on.
    
    return:
        all_tb_op: list of tuples in the form of '''
    # print(p, q, r, s)
    ordered_indices = sorted([p,q,r,s])

    tb_op_base = get_tb_base(n, ordered_indices)

    creation = [('X', 0.5), ('Y', -0.5j)]
    annihilation = [('X', 0.5), ('Y', 0.5j)]
    c_times_a = [('I', 0.5), ('Z', -0.5)]

    appeared = set()
    for i in [p,q,r,s]:    
        if i in appeared:
            repeated_index = i
        else:
            appeared.add(i)

    ops = []
    ordered_indices.remove(repeated_index) 
    for i in ordered_indices: # now len 3
        if i == repeated_index:
            ops.append(c_times_a)
        elif i == p or i == q:
            ops.append(creation)
        else:
            ops.append(annihilation)

    op_products = gen_prod_3u(ops)

    all_tb_op = []
    for op, factor in op_products:
        # every product turn into one operator in the sum
        if np.allclose(coeff, 0):
            continue # skip coeffs close to 0
        
        new_op = tb_op_base.copy()
        
        for i, index in enumerate(ordered_indices):
            new_op[index] = op[i] # op is a pauli string for the indices given not the whole qubit oprator
        
        new_op = ''.join(new_op)
        all_tb_op.append((new_op, factor * coeff))
        # print(new_op, coeff)

    return all_tb_op 



def jordan_wigner(hamiltonian):
    
    '''
    Perform the Jordan-Wigner transform for which maps second quantisation Hamiltonian to qubit Hamiltonians.
    N.B. 0th qubit to the left, always.
    arg:
        hamiltonian: Hamiltonian, the second quantisation Hamiltonian containing the overlap integrals.
    
    return:
        h_pauli: list of tuples, the qubit Hamiltonian in terms of Pauli strings.'''
    h_pauli = []
    n = hamiltonian.one_body_coeffs.shape[0]

    # one body
    for i in range(n):
        for j in range(n):

            if np.allclose(hamiltonian.one_body_coeffs[i, j],0):
                continue

            # diagnol terms
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


    if _check_health_jw(h_pauli, n):
        raise ValueError("There are terms whose length is not equal to the number of qubits.")
    
    h_pauli = _simplify_pauli_terms(h_pauli)
    return h_pauli



def _check_health_jw(h_pauli, n):
    count = 0
    for term in h_pauli:
        if len(term[0]) != n:
            print(term)
            count += 1
    
    return count
    


def _simplify_pauli_terms(terms):

    '''
        Simplify the pauli sums by collecting like terms.
        
        
        return: list, pauli sum equal to the input with like terms collected and sorted.
    
    '''
    simplified_terms = {}
    
    for pauli_string, coeff in terms:
        if pauli_string in simplified_terms:
            simplified_terms[pauli_string] += coeff
        else:
            simplified_terms[pauli_string] = coeff
    
    simplified_terms_list = [(pauli_string, coeff) for pauli_string, coeff in simplified_terms.items() if coeff != 0]
    
    # will sort alphabetically
    return sorted(simplified_terms_list)

