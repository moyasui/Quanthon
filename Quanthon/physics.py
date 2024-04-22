class Hamiltonian:

    def __init__(self, ob_coeffs, tb_coeffs) -> None:
        '''
        args:
            ob_coeffs: 2d array of one body coefficients.
            
            tb_coeffs: 4d array of two body coefficients.
        '''
        self.one_body_coeffs = ob_coeffs 
        self.two_body_coeffs = tb_coeffs

    def __repr__(self) -> str:
        return f'Hamiltonian: 1b {self.one_body_coeffs.shape}, 2b {self.two_body_coeffs.shape}.'

 