import numpy as np

def is_hermitian(mat):
    return np.allclose(mat, mat.T.conj())


def is_unitary(mat):
    return np.allclose(mat @ mat.T.conj(), np.eye(len(mat)))

def count_parametrised_gates(circuit):

    for gate in circuit:
        if gate.is_parametrised:
            count += 1
    
    return count