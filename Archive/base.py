import numpy as np
from collections import Counter

# Constants
rng = np.random.default_rng()
rangle = '\u27E9'

PendingDeprecationWarning("This module is deprecated and will be removed in the next release. Please use the new_base module instead.")

class Gate:

    def __init__(self, name, matrix, n_qubits, params=None):
        self.name = name
        self.matrix = matrix
        self.n_qubits = n_qubits
        self.params = params
    
    def __repr__(self):
        return f"Gate: {self.name} \n Matrix: \n {self.matrix} \n"
    
class Qubit:

    def __init__(self) -> None:
        self.state = np.zeros(2, dtype=np.complex_)
        self.state[0] = 1
        self.I = np.eye(2)
        self.H = 1/np.sqrt(2) * np.array([[1, 1], [1, -1]])
        self.x = np.array([[0, 1], [1, 0]])
        self.y = np.array([[0, -1j], [1j, 0]])
        self.z = np.array([[1, 0], [0, -1]])
        self.s = np.array([[1, 0], [0, 1j]])
        self.n_qubit = 1
        self._get_state_dict()
        self._get_gate_history()

    def _get_gate_history(self):
        self.gate_history = {}
        for i in range(self.n_qubit):
            self.gate_history[f'{i}'] = []
        # print(self.n_qubit)
        # print(self.gate_history)
    
    def _update_gate_history(self, gate, i):
        
        if gate in ['H', 'X', 'Y', 'Z', 'Sdag']:
            self.gate_history[f'{i}'].append(gate) 
        
        elif gate.startswith('Rx') or gate.startswith('Ry'):
            self.gate_history[f'{i}'].append(gate)
 
        elif gate == 'CNOT':
            self.gate_history[f'{i[0]}'].append(gate + f"ctrl{i[1]-i[0]}") 
            self.gate_history[f'{i[1]}'].append(gate + "trgt")
        
        elif gate == 'SWAP':
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


    def __repr__(self) -> str:
        return f"Qubit(s) in state: \n {self._to_comp_basis()} \n"

    def set_state(self, state):
        assert len(state) == 2**self.n_qubit, f"Invalid state: must have length {2**self.n_qubit}"
        assert np.linalg.norm(state) == 1, "Invalid state: must be normalised"
        self.state = state

    def copy(self):
        new_qubit = Qubit()
        new_qubit.state = self.state.copy()
        new_qubit.n_qubit = self.n_qubit
        return new_qubit

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
    
    def operate(self, op, indx): 
        ''' Operates the qubit with the given operator and index. '''
        result = np.eye(2**indx)
        result = np.kron(result, op)
        for _ in range(int(self.n_qubit-indx-len(op)/2)):
            result = np.kron(result, np.eye(2))
        
        self.state = result @ self.state
    

    def hadamard(self,i):
        self.operate(self.H,i)
        self._update_gate_history('H', i)
    
    def X(self,i):
        self.operate(self.x, i)
        self._update_gate_history('X', i)

    def Y(self,i):
        self.operate(self.y, i)
        self._update_gate_history('Y', i)

    def Z(self,i):
        self.operate(self.z, i)
        self._update_gate_history('Z', i)
    
    def sdag(self,i):
        self.operate(self.s.conj(), i)
        self._update_gate_history('Sdag', i)

    def rx(self, theta, i):
        Rx = np.cos(theta/2) * self.I - 1j * np.sin(theta/2) * self.x 
        self.operate(Rx, i)
        self._update_gate_history(f'Rx_{theta}', i)

    def ry(self, phi, i):
        Ry = np.cos(phi/2) * self.I - 1j * np.sin(phi/2) * self.y
        self.operate(Ry, i)
        self._update_gate_history(f'Ry_{phi}', i)

    def prob(self):
        prob = np.abs(self.state**2)
        return prob

    def measure(self, n_shots=1):
        ''' n: number of shots 
            indexs: the index of the qubit(s) being measured '''
        
        prob = self.prob()
        allowed_outcomes = np.arange(len(self.state))
        # print(self.state)
        outcomes = np.random.choice(allowed_outcomes, p=prob, size = n_shots)
        
        self.state = np.zeros_like(self.state)
        self.state[outcomes[-1]] = 1
        counts = Counter(outcomes)
        
        outcomes_count = np.zeros((len(self.state),2)) # 2: state and count
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
                'Rx': lambda angle: '\\gate{R_x(' + f"{angle}" + ')}',
                'Ry': lambda angle: '\\gate{R_y(' + f"{angle}" + ')}',
                'Sdag': '\\gate{S^\\dagger}',
                'CNOTctrl': lambda dist: "\\ctrl{" + f"{dist}" + '}',
                'CNOTtrgt': '\\targ{}',
                'SWAP1': lambda dist: "\\swap{" + f"{dist}" + '}',
                'SWAP2': '\\targX{}',
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

        

class Qubits_2(Qubit):

    def __init__(self):
        super().__init__()
        self.state = np.zeros(4, dtype=np.complex_)
        self.state[0] = 1
        self.n_qubit = 2
        self._get_state_dict()
        self.cnot_01 = np.array([[1, 0, 0, 0], 
                                [0, 1, 0, 0], 
                                [0, 0, 0, 1], 
                                [0, 0, 1, 0]])
        
        self.cnot_10 = np.array([[1, 0, 0, 0], 
                                [0, 0, 0, 1], 
                                [0, 0, 1, 0], 
                                [0, 1, 0, 0]])
        
        self.swp = np.array([[1, 0, 0, 0], 
                              [0, 0, 1, 0], 
                              [0, 1, 0, 0], 
                              [0, 0, 0, 1]])
        self._get_gate_history()

    def copy(self):
        new_qubit = Qubits_2()
        new_qubit.state = self.state.copy()
        new_qubit.n_qubit = self.n_qubit
        return new_qubit

    def cnot(self, i, j):
        if i == 0 and j == 1:
            self.state = self.cnot_01 @ self.state
        elif i == 1 and j == 0:
            self.state = self.cnot_10 @ self.state
        else:
            raise ValueError("Invalid indices")
        
    
    def swap(self, qubit_i=0, qubit_j=1):
        self.state = self.swp @ self.state 


class Qubits(Qubits_2):

    def __init__(self,n):
        super().__init__()
        self.state = np.zeros(2**n, dtype=np.complex_)
        self.state[0] = 1
        self.n_qubit = n
        self._get_state_dict()
        # print(self.state_dict)
        self._get_gate_history()
        
    
    def copy(self):
        new_qubits = Qubits(self.n_qubit)
        new_qubits.state = self.state.copy()
        return new_qubits


    def _find_flipped_state(self, target, state):
        ''' 
        input: target: int (0,n_qubit-1), index of the target qubit. 
        Given a state, |abcd>, find the index where the target qubit is flipped. '''
        flipped_state = list(state)
        if state[target] == "1":
            flipped_state[target] = '0'
        else:
            flipped_state[target] = '1'

        return "".join(flipped_state)
    
    def cnot(self, control, target):
        
        if self.n_qubit == 1:
            raise ValueError("The CNOT gate can not be applied to a single qubit.")
        self._get_state_dict()
        new_state_dict = self.state_dict.copy()
        # print(new_state_dict)
        for state in self.state_dict.keys():
            if state[control] == "1":
                flipped_state = self._find_flipped_state(target, state)
                new_state_dict[flipped_state] = self.state_dict[state]

        # print(len(new_state_dict.values()))
        self.state = np.fromiter(new_state_dict.values(), dtype=np.complex_)
        self._update_gate_history("CNOT", (control, target))

    def _find_swapped_state(self, i, j, state):
        swapped_state = list(state)

        swapped_state[i] = state[j] 
        swapped_state[j] = state[i]

        return "".join(swapped_state)

        
    def swap(self, qubit1, qubit2):

        if self.n_qubit == 1:
            raise ValueError("The SWAP gate can not be applied to a single qubit.")
        self._get_state_dict() 
        new_state_dict = self.state_dict.copy()

        for state in self.state_dict.keys():
            swapped_state = self._find_swapped_state(qubit1, qubit2, state)
            new_state_dict[swapped_state] = self.state_dict[state]

        # print(len(new_state_dict.values()))
        self.state = np.fromiter(new_state_dict.values(), dtype=np.complex_)
        self._update_gate_history("SWAP", (qubit1, qubit2))
    


if __name__ == "__main__":
    qc = Qubits(4)
    print(qc.n_qubit)
    print(qc)

    qc.hadamard(0)
    qc.X(1)
    qc.cnot(1,0)
    # qc.swap(0,1)
    qc.Y(0)
    # qc.rx(0.4, 1)
    
    qc.rx(-1.2080928149562626, 1)
    print(qc.gate_history)
    qc.draw(True)
