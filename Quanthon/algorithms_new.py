'''Doc:
    Normal VQE, Adapt-VQE
'''

from scipy.optimize import minimize
import numpy as np
import warnings
from .expectation import cal_expectation
from .utils import pauli_sum
from .new_base import Qubits

class VQE():
        def __init__(self, ansatz, init_points, optimiser=minimize):
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

            self.minimise = optimiser
            
                
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
     
    def __init__(self, ansatz, expectation=None, optimiser=minimize, estimate_energy=True):
        '''
        args:
            ansatz: An instance of QubitAdaptAnsatz to be optimised,
            expectation: the function which estimates the energy given a Hamiltonian,
            optimiser: the minimisation function used, default to be scipy.opmtimize.minimize,
            estimate_energy: if True, uses the estimator and shots, otherwise uses the matrix product 
                             of the state and the Hamiltonian.
        '''
        self.old_params = [] # list of parameters for the adapt circuit 
        # self.old_gates = []
        self.ansatz = ansatz
        self.expectation = cal_expectation

        self.minimise = optimiser
        
    
    def gradient(self, A):
        
        '''
        A: the operator whose gradient is calculated.
        '''
        # in practice one should estimate this as well
        # print(f"params: {self.params}")
        # print(self.ansatz.qubits)
        state = self.ansatz.run_without_update(self.params)
        # print(self.ansatz.qubits)
        
        # apply all gates first

        grad = state.conj().T @ (self.H_mat @ A - A @ self.H_mat) @ state

        
        return grad

    def _append_operator(self, eps=10e-2):
        
        '''
            Append a new operator to the ansatz which has the largest greadient ~ [H, A]
            return:
                True if largest gradient is < eps and a new operator is added, False otherwise.
        '''

        op_cand = self.ansatz.pool[0]
        max_abs_grad = 0

        # print(op_cand)
        
        for op in self.ansatz.pool:
            op_mat = pauli_sum([(op.strip('i'), 1j)])
            op_grad = self.gradient(op_mat)
            abs_grad = np.abs(op_grad)
            
            print(f"op: {op}, grad: {op_grad}")

            if abs_grad > max_abs_grad:
                op_cand = op
                max_abs_grad = abs_grad
        
        if max_abs_grad < eps:
            print(f"end of adapt, grad = {max_abs_grad}")
            return False # gradient = 0, end of adapt 
        
        # print(f"appendding, op: {op_cand}, grad: {max_abs_grad}")
        # grow the ansatz
        # self.old_gates.append(op_cand) # type(op_cand) == str
        self.params = np.append(self.params, 0)
        self.ansatz.append_op(op_cand)

        return True
    
    def _objective(self, params):
        # print("called _objective")
        state = self.ansatz.run_without_update(params)
        qc = Qubits(self.ansatz.qubits.n_qubit)
        qc.set_state(state)

        energy = self.expectation(qc, self.H_str, self.num_shots)
        # energy = state.conj().T @ self.H_mat @ state
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
        self.params = np.array([])
        for i in range(max_iter): 
            # need to update the state too
            print(self.ansatz.qubits)

            old_params, energy = self.minimise_eigenvalue(num_shots) # state is not updated here
            self.params = old_params
            print(f"i: {i}, min_energy = {energy}")
            

            if not self._append_operator(): # new parameter is added here 
                break
            
        self.ansatz.run(self.params) # update the state with the final parameters
        return self.ansatz, energy
    

    def minimise_eigenvalue(self, num_shots=10000):

        self.num_shots = num_shots
        
        # print(self.params)
        result = self.minimise(self._objective, self.params, method='Powell', options= {"maxiter": 100000})
        min_params = result.x
        min_energy = result.fun

        # self.ansatz.run(min_params) # update the state using only the optimal parameters
         
        return min_params, min_energy
