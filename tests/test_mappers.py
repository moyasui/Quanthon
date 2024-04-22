import unittest

class MappersTest(unittest.TestCase):
    
    def test_pauli_dc(self):
        
        from Quanthon.mappers import pauli_decomposition
        import numpy as np
        rng = np.random.default_rng(8463)
        h = np.array([[1,1],[1,2]])
        h2 = np.array([[1,2,1,2],
                       [2,2,1,1],
                       [1,1,3,2],
                       [2,1,2,4]])

        h2 += h2.T

        h_q_op = pauli_decomposition(h)
        h2_q_op = pauli_decomposition(h2)

        print(f"h_q_op: {h_q_op}")
        print(f"h2_q_op: {h2_q_op}")

    

if __name__ == '__main__':
    unittest.main()