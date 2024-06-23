class AdaptHistory:

    '''Contains results from AdaptVQE at each step.'''
    def __init__(self) -> None:
        self.min_energies = []
        self.operators_appended = []
        self.max_grads = []

    def update_energies(self, energy):
        self.min_energies.append(energy)
    
    def update_operators(self, operator):
        self.operators_appended.append(operator)
    
    def update_grad(self, grad):
        self.max_grads.append(grad)


    

    