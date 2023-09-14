import unittest
from Quanthon import Qubit, Qubits_2, Qubits


class QuanthonTest(unittest.TestCase):
    def test_qubit(self):
        q1 = Qubit()
        print(q1)
        q1.hadamard(0)
        print(q1)
    
    def test_qubit2(self):
        q2 = Qubits_2()
        print(q2)
        q2.cnot(0,1)
        print(q2)
        q2.swap()
        print(q2)
    
    def test_qubits(self):
        q4 = Qubits(4)
        print(q4)
        q4.state[8] = 1 
        print(q4)
        q4.swap(0,2)
        print(q4)


if __name__ == '__main__':
    unittest.main()