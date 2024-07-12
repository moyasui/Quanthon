import unittest
import numpy as np

# testing features
from Quanthon import VQE, AdaptVQE, jordan_wigner, Hamiltonian, QubitAdaptAnsatz, HardwareEfficientAnsatz, pauli_sum
from Quanthon.ansatzs import UCCSDAnsatz
# test tools
from .test_cmp_qk import qiskit_hardware
from .test_utils import is_hermitian


class AlgorithmTest(unittest.TestCase):


    # def test_vqe(self):

    #     op_str =[("II", -1.052373245772859),
    #                 ("IZ", 0.39793742484318045),
    #                 ("ZI", -0.39793742484318045),
    #                 ("ZZ", -0.01128010425623538),
    #                 ("XX", 0.18093119978423156),
    #             ]
        
    #     ansatz = HardwareEfficientAnsatz(2, reps=2)

    #     rng = np.random.default_rng(826)
    #     n_params = ansatz.n_params

    #     init_points = rng.random(n_params) * 2 * np.pi - np.pi 
    #     vqe = VQE(ansatz, init_points)
    #     min_params, min_energy = vqe.minimise_eigenvalue(op_str, 1000)
    #     print(min_energy)

    #     # with qiskit

    #     result = qiskit_hardware(op_str)
    #     # print("diff from qiskit")
    #     print(result.eigenvalue)
    
    # def _test_adapt_vqe(self):

    #     qubit_op =[("II", -1.052373245772859),
    #                 ("IZ", 0.39793742484318045),
    #                 ("ZI", -0.39793742484318045),
    #                 ("ZZ", -0.01128010425623538),
    #                 ("XX", 0.18093119978423156),
    #             ]
    #     h_mat = pauli_sum(qubit_op)
    #     # print(f"h_mat: {h_mat}")

    #     ansatz = QubitAdaptAnsatz(2)
    #     print(ansatz.qubits)

    #     print(ansatz)
        

    #     vqe = AdaptVQE(ansatz)
    #     ansatz, min_energy = vqe.run_adapt_circuit(qubit_op, 10000, max_iter=10)
    #     print(min_energy)

    #     # with qiskit

    #     result = qiskit_hardware(qubit_op)
    #     print(result.eigenvalue)

    # def _test_pairing(self):
            
    #     qubit_op = [('IIIIII', 0.1875), ('IIIIIZ', -0.5625), 
    #                 ('IIIIZI', -0.0625), ('IIIZII', 0.4375), 
    #                 ('IIZIII', -0.5625), ('IIZIIZ', 0.0625), 
    #                 ('IXXIXX', 0.03125), ('IXXIYY', (0.03125+0j)), 
    #                 ('IXYIXY', (-0.03125+0j)), ('IXYIYX', (0.03125+0j)), 
    #                 ('IYXIXY', (0.03125+0j)), ('IYXIYX', (-0.03125+0j)), 
    #                 ('IYYIXX', (0.03125+0j)), ('IYYIYY', (0.03125+0j)), 
    #                 ('IZIIII', -0.0625), ('IZIIZI', 0.0625), 
    #                 ('XXIXXI', 0.03125), ('XXIYYI', (0.03125+0j)), 
    #                 ('XYIXYI', (-0.03125+0j)), ('XYIYXI', (0.03125+0j)), 
    #                 ('XZXXZX', 0.03125), ('XZXYZY', (0.03125+0j)), 
    #                 ('XZYXZY', (-0.03125+0j)), ('XZYYZX', (0.03125+0j)), 
    #                 ('YXIXYI', (0.03125+0j)), ('YXIYXI', (-0.03125+0j)), 
    #                 ('YYIXXI', (0.03125+0j)), ('YYIYYI', (0.03125+0j)), 
    #                 ('YZXXZY', (0.03125+0j)), ('YZXYZX', (-0.03125+0j)), 
    #                 ('YZYXZX', (0.03125+0j)), ('YZYYZY', (0.03125+0j)), 
    #                 ('ZIIIII', 0.4375), ('ZIIZII', 0.0625)]
        

    #     h_mat = pauli_sum(qubit_op)
    #     if not is_hermitian(h_mat):
    #         raise Warning('This Hamiltonian is not hermitian.')

    #     my_result = self.my_hardware(qubit_op)
    #     print(my_result)

    #     result = qiskit_hardware(qubit_op)
    #     print(result.eigenvalue)

    # def my_hardware(self, qubit_op, seed=267):
    #     n = len(qubit_op[0][0])
    #     ansatz = HardwareEfficientAnsatz(n, reps=1)
    #     rng = np.random.default_rng(seed)
    #     n_params = ansatz.n_params
    #     init_points = rng.random(n_params) * 2 * np.pi - np.pi
    #     # print(init_points)
    #     vqe = VQE(ansatz, init_points)
    #     min_params, min_energy = vqe.minimise_eigenvalue(qubit_op, 10000)
    #     print(min_energy.real)


    def test_uccsd(self):

        a = UCCSDAnsatz(4, 3, 2, 2)
        a.get_singles()

if __name__ == '__main__':
    unittest.main()
