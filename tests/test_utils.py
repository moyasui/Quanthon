import numpy as np

def is_hermitian(mat):
    # print(mat == mat.T.conj())
    # for i in range(len(mat)):
    #     for j in range(len(mat)):
    #         print(i, j, mat[i][j], mat.conj().T[i][j])
    return np.allclose(mat, mat.T.conj())


def is_unitary(mat):
    return np.allclose(mat @ mat.T.conj(), np.eye(len(mat)))




