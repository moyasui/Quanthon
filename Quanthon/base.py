import numpy as np
from collections import Counter

# Constants
rng = np.random.default_rng()
rangle = '\u27E9'


class Qubit:

    def __init__(self) -> None:
        self.state = np.zeros(2, dtype=np.complex_)
        self.state[0] = 1
        self.I = np.eye(2)
        self.H = 1/np.sqrt(2) * np.array([[1, 1], [1, -1]])
        self.X = np.array([[0, 1], [1, 0]])
        self.Y = np.array([[0, -1j], [1j, 0]])
        self.Z = np.array([[1, 0], [0, -1]])
        self.S = np.array([[1, 0], [0, 1j]])
        self.n_qubit = 1
        self._get_state_dict()

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
    
    def X(self,i):
        self.operate(self.X, i)

    # Define the Pauli Y operator
    def Y(self,i):
        self.operate(self.Y, i)

    # Define the Pauli Z operator
    def Z(self,i):
        self.operate(self.Z, i)
    
    def sdag(self,i):
        self.operate(self.S.conj(), i)

    def rx(self, theta, i):
        Rx = np.cos(theta/2) * self.I - 1j * np.sin(theta/2) * self.X 
        self.operate(Rx, i)

    def ry(self, phi, i):
        Ry = np.cos(phi/2) * self.I - 1j * np.sin(phi/2) * self.Y
        self.operate(Ry, i)

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

    def copy(self):
        new_qubit = Qubits_2()
        new_qubit.state = self.state.copy()
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
        
    
    def copy(self):
        new_qubit = Qubits(self.n_qubit)
        new_qubit.state = self.state.copy()
        return new_qubit


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

        self._get_state_dict()
        new_state_dict = self.state_dict.copy()
        # print(new_state_dict)
        for state in self.state_dict.keys():
            if state[control] == "1":
                flipped_state = self._find_flipped_state(target, state)
                new_state_dict[flipped_state] = self.state_dict[state]

        # print(len(new_state_dict.values()))
        self.state = np.fromiter(new_state_dict.values(), dtype=np.complex_)

    def _find_swapped_state(self, i, j, state):
        swapped_state = list(state)

        swapped_state[i] = state[j] 
        swapped_state[j] = state[i]

        return "".join(swapped_state)

        
    def swap(self, qubit1, qubit2):

        self._get_state_dict() 
        new_state_dict = self.state_dict.copy()

        for state in self.state_dict.keys():
            swapped_state = self._find_swapped_state(qubit1, qubit2, state)
            new_state_dict[swapped_state] = self.state_dict[state]

        # print(len(new_state_dict.values()))
        self.state = np.fromiter(new_state_dict.values(), dtype=np.complex_)
        
    
    