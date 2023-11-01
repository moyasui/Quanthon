from itertools import product
import numpy as np
from Quanthon import Qubits

def get_all_pauli(n):

    '''Returns a list of all possible Pauli strings of length n (number of qubits)'''
    pauli = 'IXYZ'
    # all_paulis = list(product(pauli, pauli, repeat=int(n/2)))
    perms = [''.join(p) for p in product(pauli, repeat=n)]
    return perms


def _no_pauli_after_i(pauli_str):

    detected_I = False
    for op in pauli_str:
        if not detected_I:
            if op == 'I':
                detected_I = True
                continue
            else:
                continue
        if op != 'I':
            return False
    return True

def get_swaps(pauli_str, swap_list):

    if _no_pauli_after_i(pauli_str):
        return ''.join(pauli_str), swap_list

    found_I = False
    pauli_str = list(pauli_str)
    for i, op in enumerate(pauli_str):
        if op == 'I':
            found_I = True
            I_indx = i
            continue
        if found_I:
            if op != 'I':
                non_i_indx = i
            pauli_str[I_indx] = op 
            pauli_str[non_i_indx]= 'I'
            swap_list.append((I_indx, non_i_indx))
            found_I = False

    return get_swaps(pauli_str, swap_list)

def get_cnots(pauli_str): 

    if not _no_pauli_after_i:
        raise ValueError('Pauli string must not have any Pauli operator after I.')
    cnot_pairs = []

    for i in range(len(pauli_str)-1):
        if pauli_str[i+1] == 'I':
            break
        cnot_pairs.append((i+1, i))
        
    return cnot_pairs

def change_basis(pauli_str):

    for i in pauli_str:
        if i not in 'IXYZ':
            raise ValueError(f'Invalid Pauli operator: {i}')

    swaped_pauli, swaps = get_swaps(pauli_str, [])
    cnot_pairs = get_cnots(swaped_pauli)

    change_basis_op = []

    for pair in swaps:
        change_basis_op.append(('swap', pair))
    

    for i, op in enumerate(pauli_str):
        if op == 'X':
            change_basis_op.append(('H', i))
        elif op == 'Y':
            change_basis_op.append(('Sdag', i))
            change_basis_op.append(('H', i))
        else:
            change_basis_op.append(('I', i))
 
    
    for pair in cnot_pairs:
        change_basis_op.append(('cnot', pair))
    
    return change_basis_op

def rotate_basis(qc, ops):
    '''rotate to the measurement basis'''
    
    for op_type, idx in ops:
        if op_type == 'cnot':
            qc.cnot(idx[0], idx[1])
        elif op_type == 'swap':
            qc.swap(idx[0], idx[1])
        elif op_type == 'X':
            qc.X(idx)
        elif op_type == 'Y':
            qc.Y(idx)
        elif op_type == 'Z':
            qc.Z(idx)
        elif op_type == 'H':
            qc.hadamard(idx)
        elif op_type == 'Sdag':
            qc.sdag(idx)
        elif op_type == 'I':
            continue
        else:
            raise ValueError(f'Invalid operator: {op_type}')
    

def find_state_eigval(pauli_str):

    eigenvalues = {
        'I': [1, 1],
        'X': [1, -1],
        'Y': [1, -1],
        'Z': [1, -1]
    }

    n_qubits = len(pauli_str)
    
    # Calculate total number of states
    num_states = 2**n_qubits
    state_eigenvalues = []

    for state_num in range(num_states):
        state_str = format(state_num, f'0{n_qubits}b')  # Convert to binary string
        eigenvalue = 1  # Initialize eigenvalue
        
        # Iterate over each qubit in the state
        for i in range(n_qubits):
            pauli_op = pauli_str[i]  # Corresponding Pauli operator
            qubit_val = int(state_str[i])  # Value of the qubit (0 or 1)
            eigenvalue *= eigenvalues[pauli_op][qubit_val]  # Update eigenvalue for this qubit
        
        state_eigenvalues.append(eigenvalue)
    
    return np.array(state_eigenvalues)
    
    
def expectation(qc, pauli_ops, n_shots=10000):
    
    '''
    return:
        the expectation value
    '''

    expectation = 0
    for pauli_str, coeff in pauli_ops:

        print(pauli_str,coeff)
        # state_eigval = find_state_eigval(pauli_str[::-1])
        # state_eigval = find_state_eigval(pauli_str)
        state_eigval = np.ones(len(qc.state))
        state_eigval[int(0.5*len(qc.state)):] *= -1
        print(state_eigval)
        qc_copy = qc.copy()
        cb_ops = change_basis(pauli_str)
        print(cb_ops)

        rotate_basis(qc_copy, cb_ops)

        counts = qc_copy.measure(n_shots)[:, 0]
        print("this count", counts)
        expectation += coeff * np.sum(counts * state_eigval) / n_shots
    
    # print(cb_ops)
    return expectation


if __name__ == '__main__':
    # Test Pauli operators
    
    # print(get_all_pauli(2))
    # pauli_str = 'YY'
    # for op in pauli_str:
    #     rotation_ops = get_rotations(op)
    #     print(rotation_ops)

    a = 'X'
    b = 'Y'
    c = 'IZ'
    d = 'ZXZZ'
    # test no_op_after_i
    
    # print(no_pauli_after_i(a), no_pauli_after_i(b), no_pauli_after_i(c))

    # print(no_pauli_after_i(d))

    # PASSED

    # test get_swaps
    
    # for i in [a, b, c, d]:
    #     result, swaps = get_swaps(i, [])
    #     result = "".join(result)
    #     print(result, swaps)
    
    # PASSED

    # test get_cnots

    # for paulis in [a, b, c, d]:
    #     result, swaps = get_swaps(paulis, [])
    #     result = "".join(result) 
    #     print(result)
    #     cnot_pairs = get_cnots(result)
    #     print(f"cnots: {paulis}", cnot_pairs)

    # PASSED

    # test change_basis

    # for i in [a, b, c, d]:
    #     big_op = change_basis(i)
    #     print(f"big_op for {i}", big_op)

    # PASSED

    # test expectation

    # a = 'X'
    # b = 'IYYI'
    # c = 'ZIII'
    # d = 'XIXI'

    b = 'XX'
    c = 'YY'
    d = 'IZ'

    # testing drawing, seems ok

    qc = Qubits(2)
    hamiltonian = [b]
    # coeffs = {b:0.5 , c:0.5 , d:0.5}

    # print(f"hamiltonian: {hamiltonian}")

    for term in hamiltonian:
        
        qc_copy = qc.copy()
        ops = change_basis(term)
        rotate_basis(qc_copy, ops)
        print(ops)
        print(qc_copy.gate_history)
        qc_copy.draw()

    # testing eigenvalues