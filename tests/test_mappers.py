import unittest

class MappersTest(unittest.TestCase):
    
    def test_pd(self):
        
        # from Quanthon.mappers import pauli_decomposition
        # from Quanthon.algorithms import VQE
        import Quanthon as qt
        import numpy as np
        rng = np.random.default_rng()
        # h = np.array([[1,1],[1,2]])
        h2 = np.array([[1,2,1,2],
                       [2,2,1,1],
                       [1,1,3,2],
                       [2,1,2,4]])

        h2 += h2.T
        eigs = np.zeros(10)
        vqe_eigs = np.zeros(10)
        n = 2
        
        for i in range(10):
            h = np.abs(rng.random((6, 6)))
            h += h.T
            eig_val, _ = np.linalg.eigh(h)
            print(f"eig_val: {min(eig_val)}")
            eigs[i] = min(eig_val)
            qubit_ops = qt.pauli_decomposition(h)
            ansatz = qt.HardwareEfficientAnsatz(3)
            vqe = qt.VQE(ansatz, rng.random(ansatz.n_params) * 2 * np.pi - np.pi)
            min_eig = vqe.minimise_eigenvalue(qubit_ops)[1]
            vqe_eigs[i] = min_eig
            print(min_eig)

        import matplotlib.pyplot as plt
        plt.plot(eigs, label='exact')
        plt.plot(vqe_eigs, label='vqe')
        plt.legend()
        plt.show()
        # h_q_op = qt.pauli_decomposition(h)
        # h2_q_op = qt.pauli_decomposition(h2)

        # print(f"h_q_op: {h_q_op}")
        # print(f"h2_q_op: {h2_q_op}")


    

if __name__ == '__main__':
    unittest.main()