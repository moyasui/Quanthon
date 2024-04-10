'''Doc:
    VQEs
'''

from scipy.optimize import minimize
import numpy as np
import warnings

class VQE():
        def __init__(self, ansatz, init_points, expectation=None, optimiser=None):
            '''ansatz: a parametriced circuit that takes parmas: theta and phi
                init_points: a list of initial points, must match the number of the parameters in the ansatz
                expectation: the function that calculates the expectation value of the ansatz'''
            self.ansatz = ansatz
            self.init_points = init_points # has to match the number of the parameters in the ansatz
            self.expectation = expectation
            try:
                ansatz(init_points)
            except ValueError:
                raise ValueError(f'The initial points ({init_points}) do not match the dimension of the ansatz.')

            if optimiser is None:
                self.minimize = minimize
            
                
        def _objective(self,params):
            qc = self.ansatz(params)

            # if self.expectation is None:
            #     energy = qc.state.conj() @ (self.H @ qc.state) # removed

            energy = self.expectation(qc, self.H, self.num_shots)
            return energy

        def minimise_eigenvalue(self, H_pauli_str, num_shots=10000):
            '''
            Rotates the parametrised circuit to find the minimised energy using classical 
            minimisation algorithms.
            Inputs:
            hamiltonian: a parametrised circuit that takes theta and phi, do not depend on lambda,
            num_shots: (int) number of shots,
            return: (float) minimised energy eigenvalues.'''
            self.H = H_pauli_str

            self.num_shots = num_shots
            result = self.minimize(self._objective, self.init_points, method='Powell', options= {"maxiter": 10000})
            min_params = result.x
            min_energy = result.fun

            return min_params, min_energy
    

    
class AdaptVQE(VQE):
     
    def __init__(self, ansatz, init_points, operator_pool=None, expectation=None, optimiser=None):
        super().__init__(ansatz, init_points, expectation, optimiser)
        if operator_pool is None:
            self.pool = operator_pool # need some kind of checks here
        
        # A complete pool requires at least 2n-1 (?) operators
        if len(self.pool) <= 2*ansatz.n_qubit - 1:
            warnings.warn("Pool incomplete, energy might not converge.")


    def _default_complete_pauli_pool(self):
        pool = set()
        n_qubit = self.ansatz.n_qubit
        for _ in range(n_qubit):
            # add all V's operators
            pass

        self.pool = pool
    

    def _append_operator(self):
        
        '''
            Append a new operator to the ansatz which has the largest greadient ~ [H, A]
        '''

        op_cand = self.pool[0]
        op_grad = 0

        for op in self.pool:
            if np.abs(self.gradient(self.H,op)) > op_grad:
                op_cand = op
                op_grad = self.gradient(self.H,op)
            
        
        new_ansatz = Ansatz(self.ansatz.gates, self.ansatz.parameters)