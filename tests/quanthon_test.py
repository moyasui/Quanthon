import unittest
# from Quanthon import Qubit, Qubits_2, Qubits
import Quanthon as qt
from Quanthon import Hamiltonian, jordan_wigner
import numpy as np




class QuanthonTest(unittest.TestCase):

    def test_Hamiltonian(self):

        print("Testing `Hamiltonian`")

        h_ij = np.array([[0, 1.0], [1.0, 0]])
        h_ijkl = np.random.rand(2, 2, 2, 2)
        # h_ijkl[0, 1, 0, 1] = 0.5  
        test = Hamiltonian(h_ij, h_ijkl)

        print('Testing `jordan_wigner`')
        test_jw = jordan_wigner(test)
        print(test_jw)
        # n = 4
        # test = Hamiltonian(np.ones((n,n)), np.ones((n,n,n,n)))
        # print(jw_test)

        print(100*'-')
    
    def test_qubits(self):

        q2 = qt.Qubits(2)
        # print(q2)
        q2.H(0)
        q2.CNOT(0,1)
        q2.run()
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
        
        print(100*'-')
    def test_CNOT(self):

        print("Testing `CNOT`")
        qc = qt.Qubits(4)

        qc.H(0)
        qc.H(1)
        qc.H(2)

        # print('CNOT01')
        qc.CNOT(0,1)
        # print('CNOT12')
        qc.CNOT(1,2)
        # print('CNOT23')
        qc.CNOT(2,3)

        qc.run()
        # print(sum(qc.state**2))


        qc.CNOT(2, 1)

        qc.run()


        prob = sum(qc.state**2)
        if not np.allclose(prob, 1, rtol=1e-05, atol=1e-08):
            raise Exception("CNOT makes weird things happen.")

        # print(qc)
        print(100*'-')


    
    def test_expectation(self):

        print("Testing `cal_expectation`")
        qc = qt.Qubits(4)

        qc.H(0)
        # print(qc)
        qc.run()
    
        # Ham = [
        # ('ZIZI', (1)),
        # ('ZZXX', (0.5)), 
        # ('YYII', (1)),
        # ('IXIX', (0.5))
        # ]
        Ham = [('IIII', 1)]

        H_mat = qt.pauli_sum(Ham)
        # print(H_mat)
        energy = qt.cal_expectation(qc, Ham, 100000)

        exact_energy = qc.state.conj().T @ H_mat @ qc.state

        print("exact1: ", exact_energy, "expectation1: ", energy)
        # if not np.allclose(exact_energy,energy,rtol=1e-02, atol=1e-03):
        #     print("exact1: ", exact_energy, "expectation1: ", energy)
        #     raise Exception("cal_expectation doesn't work.")
        print(100*'-')
    
    def test_exp_2(self):

        print("Test exp_2")
        qc = qt.Qubits(4)
        qc.H(0)
        qc.CNOT(0,1)
        qc.CNOT(1,2)
        qc.CNOT(2,3)
        qc.run()
        print(qc)
        Ham = [
            ('ZIZI', (1)),
            ('IXIY', (0.5)), 
            ('YYII', (1)),
            ('IXIX', (0.5))
            ]

        H_mat = qt.pauli_sum(Ham)
        # print(H_mat)
        ee = qt.cal_expectation(qc, Ham, n_shots=100000)

        new_state = H_mat @ qc.state
        re = qc.state.conj().T @ new_state

        print("exact2: ", re, "expectation2: ", ee)
        # if not np.allclose(re,ee,rtol=1e-02, atol=1e-03):
        #     print("exact2: ", re, "expectation2: ", ee)
        #     raise Exception("cal_expectation doesn't work.")

        print(100*'-')


    def test_exp_h2(self):

        Ham = [('IIII', (-0.8126179630230781+0j)),
        ('IIIZ', (0.17119774903433002+0j)),
        ('IIZI', (-0.22278593040418412+0j)),
        ('IIZZ', (0.120544822053018+0j)),
        ('IZII', (0.17119774903433002+0j)),
        ('IZIZ', (0.16862219158920946+0j)),
        ('IZZI', (0.165867024105892+0j)),
        ('XXXX', (0.04532220205287399+0j)),
        ('XXYY', (0.04532220205287399+0j)),
        ('YYXX', (0.04532220205287399+0j)),
        ('YYYY', (0.04532220205287399+0j)),
        ('ZIII', (-0.2227859304041841+0j)),
        ('ZIIZ', (0.165867024105892+0j)),
        ('ZIZI', (0.1743484418557565+0j)),
        ('ZZII', (0.120544822053018+0j))]

        qc = qt.Qubits(4)
        qc.H(0)
        qc.CNOT(0,1)
        qc.CNOT(1,2)
        qc.CNOT(2,3)
        qc.run()
        print(qc)
        
        H_mat = qt.pauli_sum(Ham)
        # print(H_mat)
        ee = qt.cal_expectation(qc, Ham, n_shots=100000)

        new_state = H_mat @ qc.state
        re = qc.state.conj().T @ new_state

        print("exact3: ", re, "expectation3: ", ee)

    def test_exponential(self):
        
        from Quanthon import Qubits, exponential_pauli
        from Quanthon.base import Gate
        from Quanthon.utils import get_pauli
        from scipy.linalg import expm

        n = 2
        
        pauli_str = 'YX'

        a = 0.5
        coeff = -1j * a

        # with staircase algorithm
        qc = Qubits(n)
        exponential_pauli(qc, pauli_str, a, method='staircase')
        exponential_pauli(qc, pauli_str, a, method='staircase')
        exponential_pauli(qc, pauli_str, a, method='staircase')
        qc.run()
        print("staircase", qc)

        # with inverted staircase algorithm
        qc = Qubits(n)
        exponential_pauli(qc, pauli_str, a, method='inverted staircase')
        exponential_pauli(qc, pauli_str, a, method='inverted staircase')
        exponential_pauli(qc, pauli_str, a, method='inverted staircase')
        qc.run()
        print("inverted", qc)

        # with scipy.linalg.expm
        qc = Qubits(n)
        qc.reset_circuit()
        qc.circuit.append(Gate(f'exp({pauli_str})', expm(coeff * get_pauli(pauli_str)), n_qubits=qc.n_qubit))
        qc.circuit.append(Gate(f'exp({pauli_str})', expm(coeff * get_pauli(pauli_str)), n_qubits=qc.n_qubit))
        qc.circuit.append(Gate(f'exp({pauli_str})', expm(coeff * get_pauli(pauli_str)), n_qubits=qc.n_qubit))
        qc.run()
        print("scipy", qc)

if __name__ == '__main__':
    unittest.main()