'''The basic classes for quantum computing. Replaced base.py, which was discontinued after Quanthon 0.3.0'''

import numpy as np
from collections import Counter
from .utils import one_fixed_bit, flip_bit, swap_bits, get_bit, make_op_mat, is_valid_state

# Constants
rng = np.random.default_rng()
rangle = '\u27E9'

class Gate:

    def __init__(self, name, matrix, n_qubits, is_parametrised=False):
        '''
        args:
            name: string, the name of the gate;
            matrix: operator matrix or a function which takes params is an argument and returns the matrix of the correct size;
            n_qubits: int, the number of qubits the gate acts on;
            targets: the target qubit or qubits the gate acts on.
            '''
        self.name = name
        self.matrix = matrix

        self.n_qubits = n_qubits
        self.is_parametrised = is_parametrised
        # self.params = params


    def __repr__(self):
        return f"Gate: {self.name} \n Matrix: \n {self.matrix} \n"
    
    def _check_is_unitary(self, matrix):
        
        uudag = matrix @ matrix.conj().T
        if not np.allclose(uudag, np.eye(matrix.shape[0]), atol=1e-4):
            raise ValueError(f"{self.name} is not unitary, matrix UU+: {uudag}.")


    def act(self, state, param=None):

        if self.is_parametrised:
            assert param is not None
            self._check_is_unitary(self.matrix(param))
            result = self.matrix(param) @ state
        
        else:
            self._check_is_unitary(self.matrix)
            result = self.matrix @ state
        

        if is_valid_state(result):
            result /= np.linalg.norm(result)
            return result
        else:
            raise ValueError(f"Invalid state. {result} {sum(np.abs(result**2))}")

    
class Qubits:

    def __init__(self,n):
        self.state = np.zeros(2**n, dtype=np.complex128)
        self.state[0] = 1
        self.n_qubit = n
        self.operator_size = 2**n

        self.I = np.eye(2)
        self.h = 1/np.sqrt(2) * np.array([[1, 1], [1, -1]], dtype=np.complex128)
        self.x = np.array([[0, 1], [1, 0]],dtype=np.complex128)
        self.y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
        self.z = np.array([[1, 0], [0, -1]], dtype=np.complex128)
        self.s = np.array([[1, 0], [0, 1j]], dtype=np.complex128)

        self.circuit = []

        self._get_state_dict()
        self._get_gate_history()
        
    def __repr__(self) -> str:
        return f"Qubit(s) in state: \n {self._to_comp_basis()} \n"

    def _get_gate_history(self):
        self.gate_history = {}
        for i in range(self.n_qubit):
            self.gate_history[f'{i}'] = []
        # print(self.n_qubit)
        # print(self.gate_history)
    
    def _update_gate_history(self, gate, i):
        
        if gate in ['H', 'X', 'Y', 'Z', 'Sdag']:
            self.gate_history[f'{i}'].append(gate) 
        
        elif gate.startswith('R'):
            self.gate_history[f'{i}'].append(gate)
 
        elif gate.startswith('CNOT'):
            self.gate_history[f'{i[0]}'].append(gate + f"ctrl{i[1]-i[0]}") 
            self.gate_history[f'{i[1]}'].append(gate + "trgt")
        
        elif gate.startswith('SWAP'):
            self.gate_history[f'{i[0]}'].append(gate + f"1{i[1]-i[0]}")
            self.gate_history[f'{i[1]}'].append(gate + "2")

        else:
            raise ValueError(f"Invalid gate {gate}")
        
        for j in range(self.n_qubit):
            if isinstance(i, tuple):
                if j not in i:
                    self.gate_history[f'{j}'].append('I')
            else:
                if j != i:
                    self.gate_history[f'{j}'].append('I')


    def set_state(self, state):
        assert len(state) == 2**self.n_qubit, f"Invalid state: must have length {2**self.n_qubit}"
        assert np.isclose(np.linalg.norm(state), 1, atol=1e-4), "Invalid state: must be normalised."
        self.state = state

    def reset_state(self):
        self.state = np.zeros(2**self.n_qubit, dtype=np.complex_)
        self.state[0] = 1
    
    def reset_circuit(self):
        self.circuit = []

    def copy(self):
        new_qubits = Qubits(self.n_qubit)
        new_qubits.state = self.state.copy()
        return new_qubits


    def _get_state_dict(self):
        self.state_dict = dict()
        for i in np.arange(2 ** self.n_qubit):
            self.state_dict[format(i, f'0{self.n_qubit}b')] = self.state[i]
        
    def _to_comp_basis(self):
        # Calculate the number of qubits'

        # Find the indices of the states with non-zero amplitudes
        non_zero_indices = np.where(np.abs(self.state) > 0)[0]

        # Convert each non-zero index to binary representation, along with the corresponding coefficient
        computational_str = ""
        np.zeros(2**self.n_qubit, dtype=np.complex_)
        for idx in non_zero_indices:
            binary_str = format(idx, f"0{self.n_qubit}b") # set to binary
            amplitude = self.state[idx]
            computational_str += f"{amplitude:.2f}|{binary_str}⟩ + "

        computational_str = computational_str.rstrip("+ ")

        # Return the computational basis notation
        return computational_str
    
    def count_gates(self):
        return len(self.circuit)
    
    def count_cnots(self):
        return sum([1 for gate in self.circuit if gate.name.startswith('CNOT')])

    def H(self,i):
        # self.operate(self.h,i)
        matrix = make_op_mat(self.n_qubit, self.h, i)
        self.circuit.append(Gate(f'H_{i}', matrix, self.n_qubit))
        self._update_gate_history('H', i)
    
    def X(self,i):
        # self.operate(self.x, i)
        matrix = make_op_mat(self.n_qubit, self.x, i)
        self.circuit.append(Gate(f'X_{i}', matrix, self.n_qubit))
        self._update_gate_history('X', i)

    def Y(self,i):
        # self.operate(self.y, i)
        matrix = make_op_mat(self.n_qubit, self.y, i)
        self.circuit.append(Gate(f'Y_{i}', matrix, self.n_qubit))
        self._update_gate_history('Y', i)

    def Z(self,i):
        # self.operate(self.z, i)
        matrix = make_op_mat(self.n_qubit, self.z, i)
        self.circuit.append(Gate(f'Z_{i}', matrix, self.n_qubit))
        self._update_gate_history('Z', i)
    
    def Sdag(self,i):
        # self.operate(self.s.conj(), i)
        matrix = make_op_mat(self.n_qubit, self.s.conj(), i)
        self.circuit.append(Gate(f'Sdag_{i}', matrix, self.n_qubit))
        self._update_gate_history('Sdag', i)

    def Rx(self, theta, i):

        def rx(theta):
            return np.cos(theta/2) * self.I - 1j * np.sin(theta/2) * self.x
        
        if theta is None:
            self.circuit.append(Gate(f'Rx_{i}', 
                                     lambda param: make_op_mat(self.n_qubit, rx(param), i), 
                                     self.n_qubit, 
                                     is_parametrised=True))
            
            self._update_gate_history(f'Rx_θ', i)
        else:
            matrix = make_op_mat(self.n_qubit, rx(theta), i)
            self.circuit.append(Gate(f'Rx_{theta}_{i}', matrix, self.n_qubit))

            self._update_gate_history(f'Rx_{theta}', i)

    def Ry(self, theta, i):

        def ry(theta):
            return np.cos(theta/2) * self.I - 1j * np.sin(theta/2) * self.y
        
        if theta is None:
            self.circuit.append(Gate(f'Ry_{i}', 
                                     lambda param: make_op_mat(self.n_qubit, ry(param), i), 
                                     self.n_qubit, 
                                     is_parametrised=True))
            
            self._update_gate_history(f'Ry_θ', i)
        else:
            matrix = make_op_mat(self.n_qubit, ry(theta), i)
            self.circuit.append(Gate(f'Ry_{theta}_{i}', matrix, self.n_qubit))

            self._update_gate_history(f'Ry_{theta}', i)

    def Rz(self, theta, i):

        def rz(theta):
            return np.cos(theta/2) * self.I - 1j * np.sin(theta/2) * self.z
        
        if theta is None:
            self.circuit.append(Gate(f'Rz_{i}', 
                                     lambda param: make_op_mat(self.n_qubit, rz(param), i), 
                                     self.n_qubit, 
                                     is_parametrised=True))
            
            self._update_gate_history(f'Rz_θ', i)
        else:
            matrix = make_op_mat(self.n_qubit, rz(theta), i)
            self.circuit.append(Gate(f'Rz_{theta}_{i}', matrix, self.n_qubit))

            self._update_gate_history(f'Rz_{theta}', i)
 
    def CNOT(self, control, target):

        if self.n_qubit == 1:
            raise ValueError("The CNOT gate can not be applied to a single qubit.")
        
        matrix = np.eye(self.operator_size)

        indices = one_fixed_bit(self.n_qubit, self.n_qubit - control - 1, is_decimal=True) # we do n-c cuz our 0th bit is on the left

        for i in indices:
            # print(i)
            f = flip_bit(i, self.n_qubit - target - 1) # same reason
            # print(f)
            matrix[i, i] = 0
            matrix[f, i] = 1
            matrix[i, f] = 1

        self.circuit.append(Gate(f'CNOT{control}{target}', matrix, self.n_qubit))
        # self.state = matrix @ self.state
        self._update_gate_history("CNOT", (control, target))


    def SWAP(self, qubit1, qubit2):
        if self.n_qubit == 1:
            raise ValueError("The SWAP gate can not be applied to a single qubit.")
        

        matrix = np.zeros((self.operator_size, self.operator_size))
        for i in range(self.operator_size):
            j = swap_bits(i, self.n_qubit - qubit1 - 1, self.n_qubit - qubit2 - 1)
            matrix[i, j] = 1
            matrix[j, i] = 1

        # self.state = matrix @ self.state
        self.circuit.append(Gate(f'SWAP{qubit1}{qubit2}', matrix, self.n_qubit))
        self._update_gate_history("SWAP", (qubit1, qubit2))
    
    def run(self, params=None):
        
        '''
        Execute the circuit. 
        args:
            params: the parameter to be inserted into the parametrised gates, must be in order of the gates in the circuit
            
        '''
        param_indx = 0
        for gate in self.circuit:
            if gate.is_parametrised:
                self.state = gate.act(self.state, params[param_indx])
                param_indx += 1
            else:
                self.state = gate.act(self.state)
            
         
    def run_and_reset(self):
        '''Execute the circuit and reset the circuit. Return the result.'''
        self.run()
        state = self.state
        self.reset_circuit()
        return state

    def prob(self):
        prob = np.abs(self.state**2)
        
        return prob

    def measure(self, n_shots=1):
        ''' n: number of shots 
            indexs: the index of the qubit(s) being measured '''
        
        n_shots = int(n_shots)
        self.atol = 1/n_shots
        
        prob = self.prob()
        allowed_outcomes = np.arange(len(self.state))
        # print(self.state)
        try:
            outcomes = np.random.choice(allowed_outcomes, p=prob, size = n_shots)
        except ValueError:
            if np.isclose(np.sum(prob), 1, atol=self.atol):
                new_prob = prob / np.sum(prob)
                outcomes = np.random.choice(allowed_outcomes, p=new_prob, size = n_shots) 
            else:
                raise ValueError(f"Prob: {prob}, {np.sum(prob)}.")
        
        # collapse to the measurement state
        self.state = np.zeros_like(self.state)
        self.state[outcomes[-1]] = 1
        
        counts = Counter(outcomes)
        
        outcomes_count = np.zeros((len(self.state),2), dtype=object) # 2: state and count
        for i in range(len(self.state)):
            if i in counts:
                outcomes_count[i] = counts[i], format(i, f"0{self.n_qubit}b")
            else:
                outcomes_count[i] = 0,format(i, f"0{self.n_qubit}b")
        # count_qubit = Counter(outcomes_count,)
        # single_qubit_count = Counter(outcomes_count)

        return outcomes_count
        # return outcomes
    
    def draw(self, use_quantikz=True):
        print("cricuit")
        if not use_quantikz:
            raise NotImplementedError("This feature is not implemented yet.")
        else:
            gate_map = {
                'H': '\\gate{H}',
                'I': '\\qw',
                'X': '\\gate{X}',
                'Y': '\\gate{Y}',
                'Z': '\\gate{Z}',
                'Rx': lambda angle: '\\gate{R_x(' + f"{float(angle).real:.3f}" + ')}',
                'Ry': lambda angle: '\\gate{R_y(' + f"{float(angle).real:.3f}" + ')}',
                'Sdag': '\\gate{S^\\dagger}',
                'CNOTctrl': lambda dist: "\\ctrl{" + f"{dist}" + '}',
                'CNOTtrgt': '\\targ{}',
                'SWAP1': lambda dist: "\\swap{" + f"{dist}" + '}',
                'SWAP2': '\\targX{}',
                # 'exp': lambda axis, theta: f'\\gate{R_{axis}}'
            }

            max_gate_length = max(len(gates) for gates in self.gate_history.values())
            
            quantikz_str = "\\begin{quantikz}\n"
            
            for qubit in range(self.n_qubit):
                gates = self.gate_history.get(str(qubit), [])
                quantikz_str += "\t\\lstick{$q_{" + f"{qubit}" + "}$} & "
                
                for gate in gates:
                    if gate.startswith('CNOTctrl'):
                        # print(gate_map.get('CNOTctrl'))
                        quantikz_str += gate_map.get('CNOTctrl')(gate[8:])+ " & "
                    elif gate.startswith('SWAP1'):
                        quantikz_str += gate_map.get('SWAP1')(gate[5:]) + " & "
                    elif gate.startswith('Rx'):
                        quantikz_str += gate_map.get('Rx')(gate[3:]) + " & "
                    elif gate.startswith('Ry'):
                        quantikz_str += gate_map.get('Ry')(gate[3:]) + " & "
                    else:
                        quantikz_str += gate_map.get(gate, gate) + " & "
                quantikz_str += "\qw \\\\\n"
            
            # Fixing the positions of CNOT gates
            # quantikz_str = quantikz_str.replace('\\ctrl & \\targ', '\\targ & \\ctrl')
            
            quantikz_str += "\\end{quantikz}"
            
            print(quantikz_str)


if __name__ == "__main__":
    qc = Qubits(2)
    # qc.H(1)
    # qc.H(0)

    # qc.CNOT(0,1)
    # qc.CNOT(1,2)
    # qc.CNOT(2,3)


    # qc.run()
    # print(qc)

    # qc.X(1)
    # qc.Y(0)
    # qc.rx(0.4, 1)

    # print(qc.gate_history)
    # qc.draw(True)

    print(qc)