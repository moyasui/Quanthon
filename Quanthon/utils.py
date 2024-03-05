# Utils
import numpy as np

def get_pauli(basis_name):
    pauli_matrices = {'I': np.eye(2, dtype=np.complex128), 
                      'X': np.array([[0, 1], [1, 0]], dtype=np.complex128), 
                      'Y': np.array([[0, -1j], [1j, 0]]), 
                      'Z': np.array([[1, 0], [0, -1]],dtype=np.complex128)}

    # Check if basis_name is a valid combination
    valid_basis = all(ch in pauli_matrices for ch in basis_name)
    if not valid_basis:
        print(f'Invalid basis name: {basis_name}')
        return

    # Calculate the matrix for the given basis_name
    matrix = np.eye(1)
    for ch in basis_name:
        matrix = np.kron(matrix, pauli_matrices[ch])

    return matrix


def pauli_sum(lst):
    """
    Computes the Pauli sum of a list of Pauli operators. Returns the matrix.
    """
    result = 0+0j
    for item in lst:
        if len(item)!= 2:
            raise ValueError(f'must be of length 2, not {len(item)}')
        
        matrix = get_pauli(item[0])
        if matrix is None:
            raise ValueError(f'Invalid Pauli operator: {item[0]}')
        
        result += np.complex128(item[1]) * matrix 

    return result

def entropy(H, level=0):
    '''computes the von Neumann entropies for subsystem A and subsystem B
    input: 
        H: the Hamiltonian matrix (with lmb already included)
    return:
        level: int, for which energy level to compute the entropy, default is 0
        S_A, S_B, the von Neumann entropies for subsystem A and subsystem B'''
    eigenvalues, eigenvectors = np.linalg.eigh(H)

    def _compute_entropy(rho):
        eigenvalues, _ = np.linalg.eig(rho)
        entropy = -np.sum(eigenvalues * np.log2(eigenvalues + 1e-10))
        return entropy
    # index of the lowest eigenvalue and eigenvector
    
    sorted_index = np.argsort(eigenvalues)
    if type(level) == int:
        level = [level]

    entropy_A = np.zeros(len(level))
    entropy_B = np.zeros(len(level))
    for i in level:
        alpha = eigenvectors[:, sorted_index[i]]

        # density matrix for the lowest energy state
        rho_0 = np.outer(alpha, np.conj(alpha))

        # partial trace of both of the subsystems
        rho_A = np.einsum("ijkl,jl->ik", rho_0.reshape(2, 2, 2, 2), np.eye(2))
        rho_B = np.einsum("ijkl,ik->jl", rho_0.reshape(2, 2, 2, 2), np.eye(2))
        
        entropy_A[i] = _compute_entropy(rho_A)
        entropy_B[i] = _compute_entropy(rho_B)

    # Print the von Neumann entropies for subsystem A and subsystem B
    return entropy_A, entropy_B


def sum_of_pow_2(a, b):
    '''generate a sequence of sums of powers of 2 from a to b'''

    return sum([2**i for i in range(a, b+1)])

def one_fixed_bit(n,c,is_decimal=False):
    '''generate the bit string of length n such that the cth bit is always 1.

    args:
        n: int, length of the bit string.
        c: the index of the bit to be fixed to 1
        
    return:
        combinations: string, list of all the bit stirngs of length n whose cth bit is 1.
    
    N.B. 
        to convert to int, specify int(bit_sting, 2) for binary input.'''

    combinations = []
    for i in range(2 ** n):
        bit_string = f'{i:0{n}b}'
        if bit_string[n-c-1] == '1': # careful '1' is not 1
            if is_decimal:
                combinations.append(i)
            else:
                combinations.append(bit_string)

    return combinations

def flip_bit(i, t):
    '''flip the t'th bit of the binary representation of i.
    args: 
        i: int
        t: int
    
    return:
        int, the new integer after flipping the t'th bit.
    '''

    return i ^ (1 << t)


def swap_bits(val, i, j): # https://stackoverflow.com/questions/12173774/how-to-modify-bits-in-an-integer
    """
    Given an integer val, swap bits in positions i and j if they differ
    by flipping their values, i.e, select the bits to flip with a mask.
    Since v ^ 1 = 0 when v = 1 and 1 when v = 0, perform the flip using an XOR.
    """
    if (val >> i) & 1 != (val >> j) & 1:
        mask = (1 << i) | (1 << j)
        val ^= mask

    return val



    
if __name__ == "__main__":
    # Test Pauli operators
    pauli_ops = [('IZZZ', 3), ('ZXYI', 2)]
    a = pauli_sum(pauli_ops)
    
    n = 78
    t = 0
    # fn = flip_bit(n, t)
    # sn = swap_bits(4, 0, 1)
    # print(f'Flipping the {t}th bit of {n} gives {fn}')
    # print(f'Swapping the {0}th bit and the {1} bit of {4} gives {sn}')

    c = 2
    n = 8
    a = one_fixed_bit(n, c, False)
    for i in a:
        print("i:",i)
        if i[n-c-1] != '1':
            raise Exception('WRONG')
        
        i = int(i, 2)
        flipped_i = flip_bit(int(i), c) 
        print(f'{flipped_i:0{n}b}')
        if f'{flipped_i:0{n}b}'[-(c+1)] != '0':
            raise Exception('WRONG')



