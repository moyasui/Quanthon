import numpy as np

class Hamiltonian:

    def __init__(self, ob_coeffs, tb_coeffs) -> None:
        '''
        args:
            ob_coeffs: 2d array of one body coefficients.
            
            tb_coeffs: 4d array of two body coefficients.
        '''
        self.one_body_coeffs = ob_coeffs 
        self.two_body_coeffs = tb_coeffs


def jordan_wigner(hamiltonian):
    
    h_pauli_str_with_coeff = []
    n = hamiltonian.one_body_coeffs.shape[0]

    # one body
    for i in range(n):
        for j in range(n):

            if hamiltonian.one_body_coeffs[i,j] == 0:
                continue

            if i == j:
                h_pauli_str_with_coeff.append((n*'I', 0.5 * hamiltonian.one_body_coeffs[i, j]))
                h_pauli_str_with_coeff.append(((i)*'I' + 'Z' + (n-1-i)*'I', 0.5 * hamiltonian.one_body_coeffs[i, j]))


                # print(i, j)
                # print(((i)*'I' + 'Z' + (n-1-i)*'I'))
            elif i < j:
                h_pauli_str_with_coeff.append(((i)*'I' + 'Y' + (j-1-i)*'Z' + 'X' + (n-j-1)*'I', -0.25j * hamiltonian.one_body_coeffs[i, j]))
                h_pauli_str_with_coeff.append(((i)*'I' + 'X' + (j-1-i)*'Z' + 'X' + (n-j-1)*'I', 0.25 * hamiltonian.one_body_coeffs[i, j]))
                h_pauli_str_with_coeff.append(((i)*'I' + 'Y' + (j-1-i)*'Z' + 'Y' + (n-j-1)*'I', 0.25 * hamiltonian.one_body_coeffs[i, j]))
                h_pauli_str_with_coeff.append(((i)*'I' + 'X' + (j-1-i)*'Z' + 'Y' + (n-j-1)*'I', 0.25j * hamiltonian.one_body_coeffs[i, j]))

                # print(i, j)
                # print(((i)*'I' + 'Y' + (j-1-i)*'Z' + 'X' + (n-j-1)*'I'))
                # print(((i)*'I' + 'X' + (j-1-i)*'Z' + 'X' + (n-j-1)*'I'))
                # print(((i)*'I' + 'Y' + (j-1-i)*'Z' + 'Y' + (n-j-1)*'I'))
                # print(((i)*'I' + 'X' + (j-1-i)*'Z' + 'Y' + (n-j-1)*'I'))
                

    # two body
                
    # print("unhealthy terms after one body", _check_health_jw(h_pauli_str_with_coeff, n))
    # return h_pauli_str_with_coeff

    for p in range(n):
        for q in range(p+1, n):
            for r in range(n):
                for s in range(r+1, n):

                    if hamiltonian.two_body_coeffs[p, q, r, s] == 0:
                        continue
                    
                    if (p == r and q == s) or (p == s and q == r):
                        
                        # i,j,k,l = sorted([p,q,r,s])
                        # h_pauli_str_with_coeff.append((n*'I', 0.25 * hamiltonian.two_body_coeffs[p, q, r, s]))
                        # h_pauli_str_with_coeff.append(((i-1)*'I' + 'Z' + (j-1-i)*'I' + 'I' + (n-j-1)* 'I', -0.25 * hamiltonian.two_body_coeffs[p, q, r, s]))
                        # h_pauli_str_with_coeff.append(((i-1)*'I' + 'I' + (j-1-i)*'I' + 'Z' + (n-j-1)* 'I', -0.25 * hamiltonian.two_body_coeffs[p, q, r, s]))
                        # h_pauli_str_with_coeff.append(((i-1)*'I' + 'Z' + (j-1-i)*'I' + 'Z' + (n-j-1)* 'I', -0.25 * hamiltonian.two_body_coeffs[p, q, r, s]))


                        h_pauli_str_with_coeff.append((n*'I', 0.25 * hamiltonian.two_body_coeffs[p, q, r, s]))
                        h_pauli_str_with_coeff.append(((p)*'I' + 'Z' + (q-1-p)*'I' + 'I' + (n-q-1)* 'I', -0.25 * hamiltonian.two_body_coeffs[p, q, r, s]))
                        h_pauli_str_with_coeff.append(((p)*'I' + 'I' + (q-1-p)*'I' + 'Z' + (n-q-1)* 'I', -0.25 * hamiltonian.two_body_coeffs[p, q, r, s]))
                        h_pauli_str_with_coeff.append(((p)*'I' + 'Z' + (q-1-p)*'I' + 'Z' + (n-q-1)* 'I', -0.25 * hamiltonian.two_body_coeffs[p, q, r, s]))

                        # print(p, q, r, s)
                        # print((p)*'I' + 'Z' + (q-1-p)*'I' + 'I' + (n-q-1)* 'I')
                        # print((p)*'I' + 'I' + (q-1-p)*'I' + 'Z' + (n-q-1)* 'I')
                        # print((p)*'I' + 'Z' + (q-1-p)*'I' + 'Z' + (n-q-1)* 'I')
                        

                        # print(p, q, r, s)

                    elif p < q < r < s:
                        
                        i,j,k,l = p,q,r,s
                        
                        op = ((i)*'I' + 'Y' + (j-1-i)*'Z' + 'X' + (k-1-j)*'I' + 'Y' + (l-1-k)*'Z' +'X' + (n-l-1)*'I') 
                        h_pauli_str_with_coeff.append((op, -0.125 * hamiltonian.two_body_coeffs[p, q, r, s])) 

                        op = ((i)*'I' + 'Y' + (j-1-i)*'Z' + 'X' + (k-1-j)*'I' + 'Y' + (l-1-k)*'Z' +'Y' + (n-l-1)*'I')
                        h_pauli_str_with_coeff.append((op, -0.125j * hamiltonian.two_body_coeffs[p, q, r, s]))

                        op = ((i)*'I' + 'Y' + (j-1-i)*'Z' + 'X' + (k-1-j)*'I' + 'X' + (l-1-k)*'Z' +'X' + (n-l-1)*'I')
                        h_pauli_str_with_coeff.append((op, 0.125j * hamiltonian.two_body_coeffs[p, q, r, s]))

                        op = ((i)*'I' + 'Y' + (j-1-i)*'Z' + 'X' + (k-1-j)*'I' + 'X' + (l-1-k)*'Z' +'Y' + (n-l-1)*'I')
                        h_pauli_str_with_coeff.append((op, 0.125 * hamiltonian.two_body_coeffs[p, q, r, s]))

                        op = ((i)*'I' + 'Y' + (j-1-i)*'Z' + 'Y' + (k-1-j)*'I' + 'Y' + (l-1-k)*'Z' +'X' + (n-l-1)*'I')
                        h_pauli_str_with_coeff.append((op, 0.125j * hamiltonian.two_body_coeffs[p, q, r, s]))

                        op = ((i)*'I' + 'Y' + (j-1-i)*'Z' + 'Y' + (k-1-j)*'I' + 'Y' + (l-1-k)*'Z' +'Y' + (n-l-1)*'I')
                        h_pauli_str_with_coeff.append((op, -0.125 * hamiltonian.two_body_coeffs[p, q, r, s]))

                        op = ((i)*'I' + 'Y' + (j-1-i)*'Z' + 'Y' + (k-1-j)*'I' + 'X' + (l-1-k)*'Z' +'X' + (n-l-1)*'I')
                        h_pauli_str_with_coeff.append((op, 0.125 * hamiltonian.two_body_coeffs[p, q, r, s]))

                        op = ((i)*'I' + 'Y' + (j-1-i)*'Z' + 'Y' + (k-1-j)*'I' + 'X' + (l-1-k)*'Z' +'Y' + (n-l-1)*'I')
                        h_pauli_str_with_coeff.append((op, 0.125j * hamiltonian.two_body_coeffs[p, q, r, s]))

                        op = ((i)*'I' + 'X' + (j-1-i)*'Z' + 'X' + (k-1-j)*'I' + 'Y' + (l-1-k)*'Z' +'X' + (n-l-1)*'I')
                        h_pauli_str_with_coeff.append((op, -0.125j * hamiltonian.two_body_coeffs[p, q, r, s]))

                        op = ((i)*'I' + 'X' + (j-1-i)*'Z' + 'X' + (k-1-j)*'I' + 'Y' + (l-1-k)*'Z' +'Y' + (n-l-1)*'I')
                        h_pauli_str_with_coeff.append((op, 0.125 * hamiltonian.two_body_coeffs[p, q, r, s]))

                        op = ((i)*'I' + 'X' + (j-1-i)*'Z' + 'X' + (k-1-j)*'I' + 'X' + (l-1-k)*'Z' +'X' + (n-l-1)*'I')
                        h_pauli_str_with_coeff.append((op, -0.125 * hamiltonian.two_body_coeffs[p, q, r, s]))

                        op = ((i)*'I' + 'X' + (j-1-i)*'Z' + 'X' + (k-1-j)*'I' + 'X' + (l-1-k)*'Z' +'Y' + (n-l-1)*'I')
                        h_pauli_str_with_coeff.append((op, -0.125 * hamiltonian.two_body_coeffs[p, q, r, s]))

                        op = ((i)*'I' + 'X' + (j-1-i)*'Z' + 'Y' + (k-1-j)*'I' + 'Y' + (l-1-k)*'Z' +'X' + (n-l-1)*'I')
                        h_pauli_str_with_coeff.append((op, -0.125 * hamiltonian.two_body_coeffs[p, q, r, s]))

                        op = ((i)*'I' + 'X' + (j-1-i)*'Z' + 'Y' + (k-1-j)*'I' + 'Y' + (l-1-k)*'Z' +'Y' + (n-l-1)*'I')
                        h_pauli_str_with_coeff.append((op, -0.125j * hamiltonian.two_body_coeffs[p, q, r, s]))

                        op = ((i)*'I' + 'X' + (j-1-i)*'Z' + 'Y' + (k-1-j)*'I' + 'X' + (l-1-k)*'Z' +'X' + (n-l-1)*'I')
                        h_pauli_str_with_coeff.append((op, 0.125j * hamiltonian.two_body_coeffs[p, q, r, s]))

                        op = ((i)*'I' + 'X' + (j-1-i)*'Z' + 'Y' + (k-1-j)*'I' + 'X' + (l-1-k)*'Z' +'Y' + (n-l-1)*'I')
                        h_pauli_str_with_coeff.append((op, -0.125 * hamiltonian.two_body_coeffs[p, q, r, s]))

                        # if len(op) != n:
                        #     print('pqrs',p,q,r,s)
                        #     print('ijkl',i,j,k,l)
                        #     print(op)
                        
    
    if _check_health_jw(h_pauli_str_with_coeff, n):
        raise ValueError("There are terms whose length is not equal to the number of qubits.")
    
    h_pauli_str_with_coeff = _simplify_pauli_terms(h_pauli_str_with_coeff)
    return h_pauli_str_with_coeff

def _check_health_jw(h_pauli_str_with_coeff, n):
    count = 0
    for term in h_pauli_str_with_coeff:
        if len(term[0]) != n:
            count += 1
    
    return count
    


def _simplify_pauli_terms(terms):

    simplified_terms = {}
    
    for pauli_string, coeff in terms:
        if pauli_string in simplified_terms:
            simplified_terms[pauli_string] += coeff
        else:
            simplified_terms[pauli_string] = coeff
    
    simplified_terms_list = [(pauli_string, coeff) for pauli_string, coeff in simplified_terms.items()]
    
    return simplified_terms_list


if __name__ == '__main__':


    # test = Hamiltonian(np.array([[1, 2], [3, 4]]), np.random.randint(1, 10, (2, 2, 2, 2)))
    h_ij = np.array([[0, 1.0], [1.0, 0]])  # One-body terms
    h_ijkl = np.zeros((2, 2, 2, 2))  # Initialize two-body terms matrix
    h_ijkl[0, 1, 0, 1] = 0.5  
    test = Hamiltonian(h_ij, h_ijkl)
    n = 4
    test = Hamiltonian(np.ones((n,n)), np.ones((n,n,n,n)))
    jw_test = jordan_wigner(test)
    print("total terms", len(jw_test))
    # print(jw_test)
