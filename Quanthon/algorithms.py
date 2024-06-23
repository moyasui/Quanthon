'''Doc:
    Normal VQE, Adapt-VQE
'''

from scipy.optimize import minimize, Bounds
import numpy as np
import random
import warnings
from .expectation import cal_expectation
from .utils import pauli_sum
from .base import Qubits
from .history import AdaptHistory

class VQE():
        def __init__(self, ansatz, init_points, optimiser=minimize, estimate_energy=True, init_state=None
                     ):
            '''
            args:
                ansatz: a parametriced circuit that takes parmas: theta and phi
                init_points: a list of initial points, must match the number of the parameters in the ansatz
                expectation: the function that calculates the expectation value of the ansatz'''
            self.ansatz = ansatz
            self.params = init_points # has to match the number of the parameters in the ansatz
            self.expectation = cal_expectation
            self.init_state = init_state
            try:
                ansatz.create_circuit(init_points)
            except ValueError:
                raise ValueError(f'The initial points ({init_points}) do not match the dimension of the ansatz.')

            self.minimise = optimiser
            self.estimate_energy = estimate_energy
            self.obj_calls = 0

                

        def _objective(self, params):

            # print(params)
            # params = np.real(params)
            self.obj_calls += 1
            self.ansatz.create_circuit(params.real, init_state=self.init_state)
            qc = self.ansatz.qubits
            qc.run()
            
            if self.estimate_energy:
                energy = self.expectation(qc, self.H, self.num_shots)
            else:
                energy = qc.state.conj().T @ pauli_sum(self.H) @ qc.state

            # self.H_mat = pauli_sum(self.H)
            # energy = self.ansatz.qubits.state.conj().T @ self.H_mat @ self.ansatz.qubits.state
            return energy.real

        def minimise_eigenvalue(self, H_pauli_str, num_shots=10000, method='Powell', opt_log=False):
            '''
            Rotates the parametrised circuit to find the minimised energy using classical 
            minimisation algorithms.
            args:
                H_pauli_str: the Hamiltonian of the system in transformed in terms of Pauli strings,
                num_shots: (int) number of shots,

            return: (tuple) optimal parameters, minimised energy eigenvalues.'''
            self.H = H_pauli_str

            self.num_shots = num_shots
            kfs = np.ones(len(self.params), dtype=bool)
            bnds = Bounds(lb=-np.pi *  np.ones(len(self.params)), ub=np.pi * np.ones(len(self.params)), keep_feasible=kfs)

            
            result = self.minimise(self._objective, 
                                   self.params.real, 
                                   method=method,
                                   bounds=bnds, 
                                   options= {"maxiter": 10000, 'disp': opt_log}
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
        self.obj_calls = 0
 
    
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

    def _append_operator(self, eps, decompose_exp, decompose_method):
        
        '''
            Append a new operator to the ansatz which has the largest greadient ~ [H, A]
            return:
                True if largest gradient is < eps and a new operator is added, False otherwise.
        '''


        max_abs_grad = 0
        random.shuffle(self.ansatz.pool)
        for op in self.ansatz.pool:
            op_mat = pauli_sum([(op.strip('i'), 1j)])
            op_grad = self.gradient(op_mat)
            abs_grad = np.abs(op_grad)
            
            if self.print_log:
                print(f"op: {op}, grad: {op_grad}")

            if abs_grad > max_abs_grad:
                op_cand = op
                max_abs_grad = abs_grad
        
        
        if max_abs_grad < eps:
            print(f"end of adapt, grad = {max_abs_grad}")
            return False # gradient = 0, end of adapt  
        # print(f"appending, op: {op_cand}, grad: {max_abs_grad}")
        # grow the ansatz
        # self.old_gates.append(op_cand) # type(op_cand) == str
        if self.print_log:
            print("max_abs_grad: ", max_abs_grad)

        self.params = np.append(self.params, 0)
        self.ansatz.append_op(op_cand, 
                              decompose_exp=decompose_exp, 
                              decompose_method=decompose_method,
                              print_log=self.print_log)
        

        self.hist.update_grad(max_abs_grad)
        self.hist.update_operators(op_cand)
        

        return True
    
    def _objective(self, params):
        # print("called _objective")
        self.obj_calls += 1
        state = self.ansatz.run_without_update(params.real)
        qc = Qubits(self.ansatz.qubits.n_qubit)
        qc.set_state(state)

        if self.estimate_energy:
            energy = self.expectation(qc, self.H_str, self.num_shots)
        else:
            energy = state.conj().T @ self.H_mat @ state

        return energy.real
    
    def run_adapt_circuit(self, 
                          H, 
                          num_shots=10000, 
                          max_iter=100, 
                          grad_eps=1e-3, 
                          method='COBYLA', 
                          decompose_exp=False,
                          decompose_method='inverted staircase',
                          print_log=False,
                          opt_log=False):
        '''
        args:
            H: list of 2-tuple such as [('II', 0.5)]
            num_shots: number of shots,
            max_iter: maximum number of iterations,
        
        return: 
            the adapted ansatz and the energy eigenvalues.
        '''
        self.num_shots = num_shots
        self.H_mat = pauli_sum(H)
        self.H_str = H
        self.params = np.array([])
        self.print_log = print_log
        self.opt_log = opt_log
        for i in range(max_iter): 
            if self.print_log:
                print(f"i: {i}")
            # need to update the state too
            # print(self.ansatz.qubits)
            # self.ansatz.run(self.params) # ??????

            energy = self._objective(self.params)
            if not self._append_operator(eps=grad_eps, decompose_exp=decompose_exp, decompose_method=decompose_method): 
                # new parameter is added here 
                break

            old_energy = self._objective(self.params) # calculated before optimisation, even though params are NOT changed in the minimisation
            new_params, energy = self.minimise_eigenvalue(method=method) # state is not updated here
            
            self.params = new_params
            energy_with_new_param = self._objective(self.params)
            if self.print_log:
                print(f"old_E: {old_energy}")
            self.hist.update_energies(energy)

            
            # energy_after_new_op = self._objective(self.params)
            # print(f"energy after adding op: {energy_after_new_op}")
            # print(f"{energy_after_new_op} = {energy}")
            
        if i == max_iter - 1:
            print(f"Maximum number of adapt iterations ({max_iter}) reached, exit criterion ({grad_eps}) not met. ")

        self.ansatz.run(self.params) # update the state with the final parameters
        return self.ansatz, energy
    

    def minimise_eigenvalue(self, method):


        '''
            Rotates the parametrised circuit to find the minimised energy using classical 
            minimisation algorithms.
            args:
                H_pauli_str: the Hamiltonian of the system in transformed in terms of Pauli strings,
                num_shots: (int) number of shots,

            return: (tuple) optimal parameters, minimised energy eigenvalues.
        '''
        
        
        # print(self.params)

        lb = -2 * np.pi * np.ones(len(self.params))
        ub = 2* np.pi * np.ones(len(self.params))
        kfs = np.ones(len(self.params), dtype=bool)

        bounds = Bounds(lb, ub, keep_feasible=kfs)

        # old_energy = self._objective(self.params)
        # min_energy = old_energy

        result = self.minimise(self._objective, 
                            self.params, 
                            method=method, 
                            bounds=bounds, 
                            options= {"maxiter": 10000, 'disp': self.opt_log}
                            )

        min_energy = result.fun
        min_params = result.x
        if self.print_log:
            print(min_params)

        # old_energy = self._objective(self.params)
        # print(f"old energy: {old_energy.real}")
        # print(f"min_energy: {min_energy.real}")
        # print(old_energy.real < min_energy.real)
        # if old_energy.real < min_energy.real:
        #     print("Warning: New energy is higher than the old energy.")

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