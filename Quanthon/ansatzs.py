'''Doc:
    Ansatz for VQE calculations. Can grow if needs be. Parametrise the qubits object.
    parameters are taken as arguments to the circuit at every iteration'''

# TODO: add init_state to HardwareEfficientAnsatz too

# DEBUGGER

def qprint(msg):
    
    print("QubitAdaptAnsatz")
    print(msg)


import numpy as np
from scipy.linalg import expm


from .new_base import Qubits, Gate
from .utils import pauli_sum

class Ansatz:

    '''Should work for any type of ansatz that are not evolving, do not use on its own.'''
    def __init__(self, n_qubits, reps=1) -> None:

        self.n_qubits = n_qubits
        self.reps = reps
        self.n_params = 2 * n_qubits * reps

  
class HardwareEfficietnAnsatz(Ansatz):

    def __init__(self, n_qubits, reps=1) -> None:
        super().__init__(n_qubits, reps)


    def create_circuit(self, init_params):
        
        self.qubits = Qubits(self.n_qubits)
        if len(init_params) != self.n_params:
            raise ValueError(f'Initial parameters do not match the number of parameters required for the ansatz: {len(init_params)} != {self.n_params}')
        
        reshaped_params = init_params.reshape(self.reps, 2*self.n_qubits)
        
        for r in range(self.reps):
            for i in range(self.n_qubits):
                self.qubits.Rx(reshaped_params[r, i], i)
                self.qubits.Ry(reshaped_params[r, i + self.n_qubits], i)
                if i % 2 == 0:
                    self.qubits.CNOT(i, i+1)
        

        

    def run(self):
        self.qubits.run()
    


class QubitAdaptAnsatz:

    def __init__(self, n_qubits, pool='V', init_state=None) -> None:

        '''Ansatz for the Adapt-VQE calculation.'''

        self.qubits = Qubits(n_qubits)

        if init_state is None:
            # initial state not sepcified, initialise randomly
            init_state = np.random.rand(2**n_qubits)
            init_state = init_state / np.linalg.norm(init_state) 
        
        self.init_state = init_state
        self.qubits.set_state(init_state)

        if pool == 'V':
            self.pool = self.create_complete_V_pool(n_qubits)
        
        elif pool == 'G':
            self.pool = self.create_complete_G_pool()
        
        else:
            self.pool = pool # custom pool, doens't have to be complete
            
        

    def __repr__(self) -> str:
        
        return f"QubitAdaptAnsatz of {self.qubits} qubits and pool: {self.pool}."
    
    def create_complete_V_pool(self, n):
        
        if n == 2:
            return ['iYZ', 'iIY']

        prev_set = self.create_complete_V_pool(n - 1)
        new_pool = [] # set()

        for prev_op in prev_set:
            new_op = prev_op + 'Z'
            new_pool.append(new_op)


        new_pool.append('i' + (n-1) * 'I' + 'Y')
        new_pool.append('i' + (n-2) * 'I' + 'YI')


        return new_pool
        

    def create_complete_G_pool(self):
        raise NotImplementedError('Not implemented yet')


    def append_op(self, op):

        '''
        op: string, representing one of the operators in the pool
        
        '''
        # no longer need to append all the gates since they are now saved

        qprint(f"append_op op: {op}")
        op_mat = pauli_sum([(op.strip('i'), 1j)])

        def parametrised_mat(param, op_mat):
            return expm(param * op_mat)
        
        # self.params.append(0) # the corresponding parameter for the gate, initialised to 0.
        
        adapt_gate = Gate(f'exp({op})', matrix=lambda param: parametrised_mat(param, op_mat), n_qubits=self.qubits.n_qubit)
        self.qubits.circuit.append(adapt_gate)
        
    
    def run(self, params):
        for i, gate in enumerate(self.qubits.circuit):
            # print(gate.matrix(params[i]))
            self.qubits.state = gate.act(self.qubits.state, params[i])
        
        # print(self.qubits.circuit)
    
    def run_without_update(self, params):
        '''Run the circuit without updating the state.'''
        old_state = self.qubits.state
        for i, gate in enumerate(self.qubits.circuit):
            self.qubits.state = gate.act(self.qubits.state, params[i])
        
        new_state = self.qubits.state
        self.qubits.state = old_state # return to the new state while returning the result of the circuit
        return new_state



    

        



    


