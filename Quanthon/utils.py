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
    Computes the Pauli sum of a list of Pauli operators.
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



# def is_basis_combination(basis):
#     pattern = r"^[IXZY]+$"
#     match = re.match(pattern, basis)
#     return bool(match)

# def change_basis(qc, basis):

#     if not is_basis_combination(basis):
#         raise ValueError(f"Invalid basis: {basis}")
         


if __name__ == "__main__":
    # Test Pauli operators
    pauli_ops = [('IZZZ', 3), ('ZXYI', 2)]
    a = pauli_sum(pauli_ops)
    
    print(a)