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


def check_health_jw(h_pauli, n):
    count = 0
    for term in h_pauli:
        if len(term[0]) != n:
            print(term)
            count += 1
    
    return count
    

def simplify_pauli_terms(terms):

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


def is_power_of_two(n):
    return n > 0 and (2 ** (n.bit_length() - 1)) == n


def expand_h_mat(h_mat):

    '''For matrices that are not a power of 2, expand the matrix to the nearest power of 2.
    
    args:
        h_mat: 2d numpy array, the matrix to be expanded.

    return:
        new_h_mat: 2d numpy array, the expanded matrix to size of smallest power of 2 possible. 
        The original matrix is at the top left corner.
    '''

    n = h_mat.shape[0]
    pow = 1
    while True:
        if n <= 2 ** pow:
           break 
        pow += 1
    
    new_h_mat = np.zeros((2 ** pow, 2 ** pow))
    new_h_mat[:n, :n] = h_mat
    for i in range(n, 2 ** pow):
        new_h_mat[i, i] = 99
    # print(new_h_mat)
    return new_h_mat
 

