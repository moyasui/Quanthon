import unittest
# from Quanthon import Qubit, Qubits_2, Qubits
import Quanthon as qt

class QuanthonTest(unittest.TestCase):
    def test_qubit(self):
        q1 = qt.Qubit()
        # print(q1)
        q1.hadamard(0)
        # print(q1)
    
    def test_qubit2(self):
        q2 = qt.Qubits_2()
        # print(q2)
        q2.hadamard(0)
        q2.cnot(0,1)
        # print(q2)
        q2.swap(1,0)
        # print(q2)
    
    def test_qubits(self):

        print("Testing `Qubits`")
        q1 = qt.Qubits(1)
        try:
            if q1.cnot(0,1):
                raise Exception("CNOT checks not working properly.")
        except (ValueError, IndexError):
            print("CNOT safe")
        
        try:
            if q1.swap(0,1):
                raise Exception("SWAP checks not working properly.")
        except (ValueError, IndexError):
            print("SWAP safe")
        

    
    def test_expectation(self):

        print("Testing `expectation`")
        qc = qt.Qubits(4)

        qc.hadamard(0)
        # qc.cnot(0,1)
        # qc.cnot(0,2)
       

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