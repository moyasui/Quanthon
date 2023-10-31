# Under construction

class Model():

    def __init__(self, mapper='jw') -> None:
        self.hamiltonian = None # TODO: what shuold this be
        mapper_dict = {'jw': self.jordan_wigner,
                       'other': self.other_mapper,
                       'special': self.special_mapper}
        
        if type(mapper) == str and mapper in self.mapper_dict.keys():
            self.mapper = mapper
        else:
            raise ValueError(f"Unknown mapper {mapper}.")
    
        self._special_mapping_exist = False
        

    def jordan_wigner(self):
        
    
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
        

class Lipkin(Model):

    def __init__(self, J) -> None:
        super().__init__()
        self._special_mapping_exist = True
        # self.hamiltonian = self._get_hamiltonian # some pre defined hamtiltonian
        
        if J > 0 and type(J) == int:
            self.J = J # the quasi-spin 
        else:
            raise ValueError('J must be a positive integer.')
        
    def _sepcial_mapping():
        pass
        
    def _get_hamiltonian(self, J):
        pass

    def estimator(self, qc, lmb, num_shots):
        '''Expectation value of the Hamiltonian'''
        return qc.expectation(self.Hamiltonian, lmb, num_shots)

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


