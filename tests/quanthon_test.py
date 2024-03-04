import unittest
# from Quanthon import Qubit, Qubits_2, Qubits
import Quanthon as qt
import numpy as np

class QuanthonTest(unittest.TestCase):

    def test_Hamiltonian(self):

        print("Testing `Hamiltonian`")

        test_ham = qt.Hamiltonian(np.ones((2,2)), np.ones((2,2,2,2)))
        test_jw = qt.jordan_wigner(test_ham)
        print(test_jw)


    
    def test_qubits(self):

        q2 = qt.Qubits(2)
        # print(q2)
        q2.H(0)
        q2.CNOT(0,1)
        # print(q2)

        print("Testing `Qubits`")
        q1 = qt.Qubits(1)

        try:
            if q1.CNOT(0,1):
                raise Exception("CNOT checks not working properly.")
        except (ValueError, IndexError):
            print("CNOT safe")
        
        try:
            if q1.SWAP(0,1):
                raise Exception("SWAP checks not working properly.")
        except (ValueError, IndexError):
            print("SWAP safe")
        

    
    def test_expectation(self):

        print("Testing `expectation`")
        qc = qt.Qubits(4)

        qc.H(0)

       

        H = [
        ('ZIZI', (1)),
        ('IZII', (0.5)), 
        ('YYII', (1)),
        ('IXIX', (0.5))
        ]

        H_mat = qt.pauli_sum(H)
        # print(H_mat)
        energy = qt.expectation(qc, H, 100000)

        exact_energy = qc.state.T.conj() @ H_mat @ qc.state

        print("exact: ", exact_energy, "expectation: ", energy)


if __name__ == '__main__':
    unittest.main()