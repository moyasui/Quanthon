'''Doc:
    VQEs
'''

from scipy.optimize import minimize
import numpy as np
import warnings
from .Ansatz import HardwareEfficietnAnsatz, QubitAdaptAnsatz
from .Expectation import expectation
from .utils import pauli_sum

class VQE():
        def __init__(self, ansatz, init_points, optimiser=None):
            '''ansatz: a parametriced circuit that takes parmas: theta and phi
                init_points: a list of initial points, must match the number of the parameters in the ansatz
                expectation: the function that calculates the expectation value of the ansatz'''
            self.ansatz = ansatz
            self.params = init_points # has to match the number of the parameters in the ansatz
            self.expectation = expectation
            try:
                ansatz(init_points)
            except ValueError:
                raise ValueError(f'The initial points ({init_points}) do not match the dimension of the ansatz.')

            if optimiser is None:
                self.minimize = minimize
            
                
        def _objective(self,params):
            self.ansatz.create_circuit(params)
            qc = self.ansatz.qubits
            qc.run()

            energy = self.expectation(qc, self.H, self.num_shots)
            return energy

        def minimise_eigenvalue(self, H_pauli_str, num_shots=10000):
            '''
            Rotates the parametrised circuit to find the minimised energy using classical 
            minimisation algorithms.
            Inputs:
            H_pauli_str: the Hamiltonian of the system in transformed in terms of Pauli strings,
            num_shots: (int) number of shots,
            return: (float) minimised energy eigenvalues.'''
            self.H = H_pauli_str

            self.num_shots = num_shots
            result = self.minimize(self._objective, self.init_points, method='Powell', options= {"maxiter": 10000})
            min_params = result.x
            min_energy = result.fun
            
            return min_params, min_energy
    

    
class AdaptVQE(VQE):
     
    def __init__(self, ansatz, init_points, expectation=None, optimiser=None):
        super().__init__(ansatz, init_points, expectation, optimiser)
        self.old_params = [] # list of parameters for the adapt circuit 
        self.old_gates = []

    
    def gradient(self, H, A):
        
        grad = self.ansatz.qubits.state.conj() @ (H @ A - A @ H) @ self.ansatz.qubits.state # in practice one should estimate this as well
        
        return grad

    def _append_operator(self, eps=10e-6): # TODO: fix this
        
        '''
            Append a new operator to the ansatz which has the largest greadient ~ [H, A]
        '''

        op_cand = self.ansatz.pool[0]
        op_grad = 0

        op_mat = pauli_sum([op_cand.strip('i'), 1j])

        for op in self.ansatz.pool:
            if np.abs(self.gradient(self.H,op_mat)) > op_grad:
                op_cand = op
                op_grad = self.gradient(self.H,op)
        
        if op_grad < eps:
            return 1
        
        # grow the ansatz
        self.old_gates.append(op_cand)
        self.params.append(0) 

        return 0
    
    def _objective(self, params):
        return super()._objective(params)
    
    def adapt_circuit(self, H, num_shots=10000, max_iter=100):
        self.H = H
        self.params = []
        for _ in range(max_iter):
            old_params = self.minimise_eigenvalue(H, num_shots)[0]
            if self._append_operator():
                break

        return self.ansatz