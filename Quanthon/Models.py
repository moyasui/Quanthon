# Under construction
'''This will be the last thing to be implemented with built in models, Hamiltonians and special mapping if any.'''
import numpy as np
from itertools import permutations


class Model():

    def __init__(self, integral_values, mapper='jw') -> None:
        '''
        args:
            integral_values: array of (*, 2) of integral values for the model.
        '''
        self.s = 'ad,a'
        
        mapper_dict = {'jw': self.jordan_wigner,
                       'other': self.other_mapper,
                       'special': self.special_mapper}
        
        self.states = np.zeros(np.log2(self.hamiltonian.shape[0]))

        if mapper not in self.mapper_dict.keys():
            self.mapper = mapper
        else:
            raise ValueError(f"Unknown mapper {mapper}.")
    
    # def hamiltonian(self, pq_max):
    #     assert(type(pq_max) == int)
    #     range_of_pq = np.arange(0, pq_max)
        

    def adag(self, i):
        ''' applies the creation operator on state i.
        N.B. This state is the energy state of the model, not the quantum state used in quantum computers.'''
        
        if self.states[i] == 1:
            raise ValueError(f"State {i} is already occupied.")
        self.states[i] = 1


    def a(self,i):
        ''' applies the annihilation operator on state i'''

        if self.states[i] == 0:
            raise ValueError(f"Nothing to annihilate at state {i}.")
        self.states[i] = 0

    def jordan_wigner(self):

        return 
    
    def other_mapper(self):
        raise NotImplementedError('Other mappers are not implemented yet.')

        
    def _to_Pauli_strings(self, mapper='jw'):
        if self._special_mapping_exist:
            self._special_mapping()  
            
        pass

    def estimator(self):
        ''' Genrate the correct rotation to do the measurements for Hamilronians given by the pauli strings.
        
        return:
            (list) of the rotations 
        '''

        pass

class Lipkin_Model():  
    '''The Lipkin (LMG) model.''' 
    #TODO: combine into Model at some point.
    
    def __init__(self, J=1, mapper='jw') -> None:
        super().__init__(mapper)
        self.J = J
        self.hamiltonian = self.Lipkin_Hamiltonian(J)
        self._special_mapping_exist = True

    def Lipkin_Hamiltonian(v, J=1, z_coeff=0.5):
        
        '''
        returns a string of pauli operator given the J values of the lipkin model.'''
        x_coeff = -v/2
        y_coeff = v/2

        pauli_str = []
        Z_str = 'I'*(2*J-1) + 'Z'
        Z_ops = set(permutations(Z_str))
        # print(Z_ops)
        Z_ops = [''.join(p) for p in Z_ops]
        for op in Z_ops:
            # print(op)
            pauli_str.append((op, z_coeff))

        # the X's   
        X_str = 'I'*(2*J-2) + 'XX'
        X_ops = set(permutations(X_str))
        X_ops = [''.join(p) for p in X_ops]

        for op in X_ops:
            # print(op)
            pauli_str.append((op, x_coeff))

        # the Y's 

        Y_str = 'I'*(2*J-2) + 'YY'
        Y_ops = set(permutations(Y_str))
        Y_ops = [''.join(p) for p in Y_ops]

        for op in Y_ops:
            # print(op)
            pauli_str.append((op, y_coeff))
        
        return pauli_str


class Ising(Model):

    def __init__(self, mapper='jw') -> None:
        super().__init__(mapper)


class Heisenberg(Model):

    def __init__(self, mapper='jw') -> None:
        super().__init__(mapper)
    
    
class SSH(Model):

    def __init__(self, mapper='jw') -> None:
        super().__init__(mapper)


class KitaevChain(Model):

    def __init__(self, mapper='jw') -> None:
        super().__init__(mapper)



