'''Doc:
    Normal VQE, Adapt-VQE
'''

from scipy.optimize import minimize
import numpy as np
import warnings
from .expectation import cal_expectation
from .utils import pauli_sum

class VQE():
        def __init__(self, ansatz, init_points, optimiser=None):
            '''
            args:
                ansatz: a parametriced circuit that takes parmas: theta and phi
                init_points: a list of initial points, must match the number of the parameters in the ansatz
                expectation: the function that calculates the expectation value of the ansatz'''
            self.ansatz = ansatz
            self.params = init_points # has to match the number of the parameters in the ansatz
            self.expectation = cal_expectation
            try:
                ansatz.create_circuit(init_points)
            except ValueError:
                raise ValueError(f'The initial points ({init_points}) do not match the dimension of the ansatz.')

            if optimiser is None:
                self.minimize = minimize
            else:
                self.minimize = optimiser
            
                
        def _objective(self, params):
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
            result = self.minimize(self._objective, self.params, method='Powell', options= {"maxiter": 10000})
            min_params = result.x
            min_energy = result.fun
            
            return min_params, min_energy
    

    
class AdaptVQE():
     
    def __init__(self, ansatz, expectation=None, optimiser=None):
        self.old_params = [] # list of parameters for the adapt circuit 
        self.old_gates = []
        self.ansatz = ansatz
        self.expectation = cal_expectation

        if optimiser is None:
                self.minimize = minimize
        else:
            self.minimize = optimiser

    
    def gradient(self, A):
        
        '''
        A: the operator whose gradient is calculated.
        '''
        # in practice one should estimate this as well
        state = self.ansatz.qubits.state
        print(f"H: {self.H_mat}, A: {A}")
        grad = state.conj().T @ (self.H_mat @ A - A @ self.H_mat) @ state
        
        return grad

    def _append_operator(self, eps=10e-6):
        
        '''
            Append a new operator to the ansatz which has the largest greadient ~ [H, A]
            return:
                1 if largest gradient is < eps, 0 otherwise.
        '''

        op_cand = self.ansatz.pool[0]
        op_grad = 0

        # print(op_cand)
        op_mat = pauli_sum([(op_cand.strip('i'), 1j)])

        for op in self.ansatz.pool:
            if np.abs(self.gradient(op_mat)) > op_grad:
                op_cand = op
                op_grad = self.gradient(op)
        
        if op_grad < eps:
            return 0 # gradient = 0, end of adapt 
        
        # grow the ansatz
        self.old_gates.append(op_cand)
        self.params.append(0) 

        return 1 # 
    
    def _objective(self, params):
        self.ansatz.run(params)
        qc = self.ansatz.qubits

        energy = self.expectation(qc, self.H_str, self.num_shots)
        return energy
    
    def run_adapt_circuit(self, H, num_shots=10000, max_iter=100):
        '''
        args:
            H: list of 2-tuple such as [('II', 0.5)]
            num_shots: number of shots,
            max_iter: maximum number of iterations,
        
        return: 
            the adapted ansatz and the energy eigenvalues.
        '''
        self.H_mat = pauli_sum(H)
        self.H_str = H
        self.params = []
        for _ in range(max_iter):
            old_params, energy = self.minimise_eigenvalue(num_shots)
            if not self._append_operator():
                break

        return self.ansatz, energy
    

    def minimise_eigenvalue(self, num_shots=10000):

            self.num_shots = num_shots
            result = self.minimize(self._objective, self.params, method='Powell', options= {"maxiter": 10000})
            min_params = result.x
            min_energy = result.fun
            
            return min_params, min_energy
    