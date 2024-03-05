'''Doc:
    Ansatz for VQE calculations. Can grow if needs be. Parametrise the qubits object.
    parameters are taken as arguments to the circuit at every iteration'''

import numpy as np
from scipy.linalg import expm


from .new_base import Qubits, Gate
from .utils import pauli_sum



class HardwareEfficietnAnsatz:

    def __init__(self, n_qubits, reps=1) -> None:
        '''args:
            params_init: array of parameters initial values'''
        self.qubits = Qubits(n_qubits)
        self.n_qubits = n_qubits
        self.reps = reps


    def create_circuit(self, params):
        
        self.qubits = Qubits(self.n_qubits)

        if len(params) != 2*self.reps:
            raise ValueError(f'Initial parameters do not match the number of parameters required for the ansatz: {len(params)} != {2*self.rep}')
        
        params.reshape(self.reps, 2*self.qubits.n_qubit)
        for r in range(self.reps):
            for i in range(self.qubits.n_qubit):
                self.qubits.rx(params[r, i], i)
                self.qubits.ry(params[r, i + self.qubits.n_qubit], i)
                if i % 2 == 0:
                    self.qubits.cnot(i, i+1)

    def run(self):
        self.qubits.run()

    


class QubitAdaptAnsatz:

    def __init__(self, n_qubits, pool_type) -> None:

        self.qubits = Qubits(n_qubits)
        if pool_type == 'V':
            self.pool = self.create_complete_V_pool(n_qubits)
        
        elif pool_type == 'G':
            self.create_complete_G_pool()

    def create_complete_V_pool(self, n):

        if n == 2:
            return {'iYZ', 'iIY'}

        prev_set = self.create_complete_V_pool(n - 1)
        new_pool = set()

        for prev_op in prev_set:
            new_op = prev_op + 'Z'
            new_pool.add(new_op)


        new_pool.add('i' + (n-1) * 'I' + 'Y')
        new_pool.add('i' + (n-2) * 'I' + 'YI')


        return new_pool

    
    def create_complete_G_pool(self):
        raise NotImplementedError('Not implemented yet')


    def append_op(self, op, params, old_gates):

        '''
        old_gates: matrices of the gates from before
        
        '''
        # make a circuit
        self.qubits = Qubits(self.n_qubits) 
        for gate in old_gates:
            self.qubits.circuit.append(gate)

        op_mat = pauli_sum([op.strip('i'), 1j])
        adapt_gate = Gate(f'exp({op})', expm(1j * op_mat))
        self.qubits.circuit.append(adapt_gate)

    
    def run(self):
        self.qubits.run()

        



    


