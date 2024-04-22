import numpy as np

def is_hermitian(mat):
    return np.allclose(mat, mat.T.conj())


def is_unitary(mat):
    return np.allclose(mat @ mat.T.conj(), np.eye(len(mat)))




