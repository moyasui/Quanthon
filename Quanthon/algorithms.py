'''Doc:
    Normal VQE, Adapt-VQE
'''

from scipy.optimize import minimize, Bounds
import numpy as np
import warnings
from .expectation import cal_expectation
from .utils import pauli_sum
from .base import Qubits
from .history import AdaptHistory

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

            # print(params)
            params = np.real(params)
            self.ansatz.create_circuit(params)
            qc = self.ansatz.qubits
            qc.run()
            
            energy = self.expectation(qc, self.H, self.num_shots)

            # self.H_mat = pauli_sum(self.H)
            # energy = self.ansatz.qubits.state.conj().T @ self.H_mat @ self.ansatz.qubits.state
            return energy

        def minimise_eigenvalue(self, H_pauli_str, num_shots=10000, method='powell'):
            '''
            Rotates the parametrised circuit to find the minimised energy using classical 
            minimisation algorithms.
            args:
                H_pauli_str: the Hamiltonian of the system in transformed in terms of Pauli strings,
                num_shots: (int) number of shots,

            return: (tuple) optimal parameters, minimised energy eigenvalues.'''
            self.H = H_pauli_str

            self.num_shots = num_shots

            lb = -np.pi * np.ones(len(self.params))
            ub = np.pi * np.ones(len(self.params))
            bounds = Bounds(lb, ub)
            
            if method == 'powell':
                bounds = None
            
            result = self.minimise(self._objective, 
                                   self.params, 
                                   method=method,
                                #    method='Powell',
                                #    method='Nelder-Mead',
                                    #  method='COBYLA',
                                    #  method='TNC',
                                   bounds=bounds, 
                                   options= {"maxiter": 100}
                                   )
            min_params = result.x
            min_energy = result.fun
            
            # print(min_params)
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
        self.estimate_energy = estimate_energy
        self.hist = AdaptHistory()
        
    
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

        rng = np.random.default_rng(4418)
        op_cand = rng.choice(self.ansatz.pool, 1).tolist()[0]
        max_abs_grad = 0

        # print(op_cand)
        
        for op in self.ansatz.pool:
            op_mat = pauli_sum([(op.strip('i'), 1j)])
            op_grad = self.gradient(op_mat)
            abs_grad = np.abs(op_grad)
            
            # print(f"op: {op}, grad: {op_grad}")

            if abs_grad > max_abs_grad:
                op_cand = op
                max_abs_grad = abs_grad
        
        
        if max_abs_grad < eps:
            print(f"end of adapt, grad = {max_abs_grad}")
            return False # gradient = 0, end of adapt  
        # print(f"appending, op: {op_cand}, grad: {max_abs_grad}")
        # grow the ansatz
        # self.old_gates.append(op_cand) # type(op_cand) == str
        print("max_abs_grad: ", max_abs_grad)
        self.params = np.append(self.params, 0)
        self.ansatz.append_op(op_cand)

        self.hist.update_grad(max_abs_grad)
        self.hist.update_operators(op_cand)
        

        return True
    
    def _objective(self, params):
        # print("called _objective")
        state = self.ansatz.run_without_update(params)
        qc = Qubits(self.ansatz.qubits.n_qubit)
        qc.set_state(state)

        if self.estimate_energy:
            energy = self.expectation(qc, self.H_str, self.num_shots)
        else:
            energy = state.conj().T @ self.H_mat @ state

        return energy
    
    def run_adapt_circuit(self, H, num_shots=10000, max_iter=100, grad_eps=1e-4, method='COBYLA'):
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
            # print(self.ansatz.qubits)
            # self.ansatz.run(self.params) # ??????
            old_params, energy = self.minimise_eigenvalue(num_shots, method=method) # state is not updated here
            self.params = old_params
            print(f"i: {i}, min_energy = {energy}")
            self.hist.update_energies(energy)

            if not self._append_operator(eps=grad_eps): # new parameter is added here 
                break
            
        if i == max_iter - 1:
            print("Maximum adapt iteration reached, energy did not converge.")
            # warnings.warn('Adapt circuit did not converge.')
        self.ansatz.run(self.params) # update the state with the final parameters
        return self.ansatz, energy
    

    def minimise_eigenvalue(self, num_shots=10000, method='SLSQP'):


        '''
            Rotates the parametrised circuit to find the minimised energy using classical 
            minimisation algorithms.
            args:
                H_pauli_str: the Hamiltonian of the system in transformed in terms of Pauli strings,
                num_shots: (int) number of shots,

            return: (tuple) optimal parameters, minimised energy eigenvalues.
        '''
        self.num_shots = num_shots
        
        # print(self.params)

        lb = -np.pi * np.ones(len(self.params))
        ub = np.pi * np.ones(len(self.params))
        bounds = Bounds(lb, ub)

        # old_energy = self._objective(self.params)
        # min_energy = old_energy

        if self.params.size == 0:
            result = self.minimise(self._objective, 
                                0, 
                                method=method, 
                                #    bounds=bounds,
                                options= {"maxiter": 100000}
                                )
        else:
            result = self.minimise(self._objective, 
                                self.params, 
                                method=method, 
                                #    bounds=bounds,
                                options= {"maxiter": 100000}
                                )

        min_energy = result.fun
        min_params = result.x

        old_energy = self._objective(self.params)
        print(f"old energy: {old_energy}")
        print(f"min_energy: {min_energy}")
        if old_energy < min_energy:
            warnings.warn("New energy is higher than the old energy.")

        # self.ansatz.run(min_params) # update the state using only the optimal parameters
         
        return min_params, min_energy


def qft(qubits, is_inverse=False, is_swap=True):

    def nth_root_of_unity(self, n):
        return 2j * n * np.pi / (2 ** self.qubits.n_qubit)
    
    if is_inverse:
        for i in range(qubits.n_qubit):
            for j in range(i):
                qubits.cp(-np.pi / 2 ** (i - j), j, i)
            qubits.H(i)
    else:
        for i in range(qubits.n_qubit):
            qubits.H(i)
            for j in range(i+1, qubits.n_qubit):
                qubits.CP(nth_root_of_unity(j), j, i) # phi, controlled, target
    
    if is_swap:
        for i in range(qubits.n_qubit/2):
            qubits.SWAP(i, qubits.n_qubits-1) 

def qpe():
    raise NotImplementedError("Quantum Phase Estimation is not implemented yet.")