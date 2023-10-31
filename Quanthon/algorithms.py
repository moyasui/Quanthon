'''Doc:
    VQEs'''

from scipy.optimize import minimize

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

            if self.expectation is None:
                energy = qc.state.conj() @ (self.H @ qc.state) # a bit cheating
            else:
                energy = self.expectation(qc, self.lmb, self.num_shots)
            return energy

        def minimise_eigenvalue(self, hamiltonian, lmb, num_shots=10000):
            '''
            Rotates the parametrised circuit to find the minimised energy using classical 
            minimisation algorithms.
            Inputs:
            hamiltonian: a parametrised circuit that takes theta and phi, do not depend on lambda,
            num_shots: (int) number of shots,
            return: (float) minimised energy eigenvalues.'''
            self.H = hamiltonian
            self.lmb = lmb
            self.num_shots = num_shots
            result = self.minimize(self._objective, self.init_points, method='Powell', options= {"maxiter": 10000})
            min_params = result.x
            min_energy = result.fun

            return min_params, min_energy
    
