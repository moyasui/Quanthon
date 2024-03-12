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

def get_ob_op(n, i, j, coeff):
    # when the indices are different
    base = ['' if k == i or k == j else 'Z' if i < k < j else 'I' for k in range(n)]
    all_op = []
    for op, factor in [('XX', 0.25), ('XY', 0.25j), ('YX', -0.25j), ('YY', 0.25)]:
        new_op = base.copy()
        new_op[i] = op[0]
        new_op[j] = op[1]
        all_op.append((''.join(new_op), factor * coeff))

    return all_op


def generate_products(creation, annihilation):
    """
    Generate products of creation and annihilation operators.
    """
    products = []
    for c1, coeff_c1 in creation:
        for c2, coeff_c2 in creation:
            for a1, coeff_a1 in annihilation:
                for a2, coeff_a2 in annihilation:
                    operator = c1 + c2 + a1 + a2
                    coefficient = coeff_c1 * coeff_c2 * coeff_a1 * coeff_a2
                    products.append((operator, coefficient))

    return products



def get_tb_base(n,ordered_indices):

    # get the I's and the Z's for the qubits not acted on. 
    tb_op_base = list(n * 'I')

    for i in range(n):
        if i in ordered_indices:
            continue

        if ordered_indices[0] < i < ordered_indices[1]:
            tb_op_base[i] = 'Z'
        
        if ordered_indices[2] < i < ordered_indices[3]:
            tb_op_base[i] = 'Z'
        
    return tb_op_base
    
def get_tb_op(n,p,q,r,s,coeff):
    # when all the indices are different

    # print(coeff, p, q, r, s)
    ordered_indices = sorted([p,q,r,s])
    tb_op_base = get_tb_base(n, ordered_indices)

    creation = [('X', 0.25 * coeff), ('Y', -0.25j * coeff)]
    annihilation = [('X', 0.25 * coeff), ('Y', 0.25j * coeff)]

    op_products = generate_products(creation,annihilation)

    # the first and the third term has odd number of Z(s), 
    # and only the annihilation operator changes sign.
    if ordered_indices[0] == r or ordered_indices[2] == r:
        op_products = [(op, coeff * -1) for op, coeff in op_products]
    
    if ordered_indices[0] == s or ordered_indices[2] == s:
        op_products = [(op, coeff * -1) for op, coeff in op_products]
    
    all_tb_op = []
    for op, coeff in op_products:
        # every product turn into one operator in the sum
        if np.allclose(coeff, 0, rtol=1e-04, atol=1e-05):
            continue # skip coeffs close to 0
        
        new_op = tb_op_base.copy()
        
        for i, index in enumerate(ordered_indices):
            new_op[index] = op[i]
        
        new_op = ''.join(new_op)
        all_tb_op.append((new_op, coeff))
        # print(new_op, coeff)
    return all_tb_op


def get_tb_one_set(n, p, q, r, s):

    assert len({p, q, r, s}) == 3

    tb_op_base = get_tb_base 
    print(tb_op_base)


def jordan_wigner(hamiltonian):
    
    h_pauli = []
    n = hamiltonian.one_body_coeffs.shape[0]

    # one body
    for i in range(n):
        for j in range(n):

            if np.allclose(hamiltonian.one_body_coeffs[i, j],0):
                continue

            if i == j:
                # II...III - I...IZ...I

                all_i = n * 'I'
                h_pauli.append((all_i, 0.5 * hamiltonian.one_body_coeffs[i, j]))

                z_term = list(all_i)
                z_term[i] = 'Z'
                z_term = ''.join(z_term)

                h_pauli.append((z_term, -0.5 * hamiltonian.one_body_coeffs[i, j]))
                # print(z_term, hamiltonian.one_body_coeffs[i, j])

                # print(i, j)
                # print(((i)*'I' + 'Z' + (n-1-i)*'I'))

            elif i < j:
                
                ob_op = get_ob_op(n, i, j, hamiltonian.one_body_coeffs[i, j])
                print(ob_op)
                h_pauli.extend(ob_op)
                
                # h_pauli.append(((i)*'I' + 'Y' + (j-1-i)*'Z' + 'X' + (n-j-1)*'I', -0.25j * hamiltonian.one_body_coeffs[i, j]))
                # h_pauli.append(((i)*'I' + 'X' + (j-1-i)*'Z' + 'X' + (n-j-1)*'I', 0.25 * hamiltonian.one_body_coeffs[i, j]))
                # h_pauli.append(((i)*'I' + 'Y' + (j-1-i)*'Z' + 'Y' + (n-j-1)*'I', 0.25 * hamiltonian.one_body_coeffs[i, j]))
                # h_pauli.append(((i)*'I' + 'X' + (j-1-i)*'Z' + 'Y' + (n-j-1)*'I', 0.25j * hamiltonian.one_body_coeffs[i, j]))

                # print(i, j)
                # print(((i)*'I' + 'Y' + (j-1-i)*'Z' + 'X' + (n-j-1)*'I'))
                # print(((i)*'I' + 'X' + (j-1-i)*'Z' + 'X' + (n-j-1)*'I'))
                # print(((i)*'I' + 'Y' + (j-1-i)*'Z' + 'Y' + (n-j-1)*'I'))
                # print(((i)*'I' + 'X' + (j-1-i)*'Z' + 'Y' + (n-j-1)*'I'))
                

    # two body
                
    # print("unhealthy terms after one body", _check_health_jw(h_pauli, n))
    # return h_pauli

    for p in range(n):
        for q in range(n):
            if p == q:
                continue
            for r in range(n):
                for s in range(n):
                    if r == s:
                        continue

                    if np.allclose(hamiltonian.two_body_coeffs[p, q, r, s],0):
                        # print(hamiltonian.two_body_coeffs[p, q, r, s])
                        continue
                    
                    # print(f'pqrs: {p}{r}{q}{s}, coeff: {hamiltonian.two_body_coeffs[p, q, r, s]}')


                    if (p == r and q == s) or (p == s and q == r):
                        
                        # print(p, q, r, s)
                        all_i = n*'I'

                        pz = list(n*'I')
                        pz[p] = 'Z'
                        pz = ''.join(pz)

                        qz = list(n*'I')
                        qz[q] = 'Z'
                        qz = ''.join(qz)

                        qpz = list(n*'I')
                        qpz[p] = 'Z'
                        qpz[q] = 'Z'
                        qpz = ''.join(qpz)

                        # print(p, q, r, s)
                        # print(all_i, pz, qz, qpz)
                        # print("all_i", all_i)
                        # print("pz", pz)
                        # print("qz", qz)
                        # print("qpz", qpz)

                        h_pauli.append((all_i, 0.25 * (hamiltonian.two_body_coeffs[p, q, r, s])))
                        h_pauli.append((pz, -0.25 * (hamiltonian.two_body_coeffs[p, q, r, s])))
                        h_pauli.append((qz, -0.25 * (hamiltonian.two_body_coeffs[p, q, r, s])))
                        h_pauli.append((qpz, 0.25 * (hamiltonian.two_body_coeffs[p, q, r, s])))

                        
                    elif len({p, q, r, s}) == 4:
                        # continue
                        # if they are all different
                        # print('pqrs',p,q,r,s)

                        tb_op = get_tb_op(n, p, q, r, s, hamiltonian.two_body_coeffs[p, q, r, s])
                        h_pauli.extend(tb_op)

                    elif len({p, q, r, s}) == 3:
                        # 0223
                        raise NotImplementedError("Not implemented yet")
                        # get_tb_one_set(n, p, q, r, s)
                        


    if _check_health_jw(h_pauli, n):
        raise ValueError("There are terms whose length is not equal to the number of qubits.")
    
    h_pauli = _simplify_pauli_terms(h_pauli)
    return h_pauli



def _check_health_jw(h_pauli, n):
    count = 0
    for term in h_pauli:
        if len(term[0]) != n:
            print(term)
            count += 1
    
    return count
    


def _simplify_pauli_terms(terms):

    '''
        Simplify the pauli sums by collecting like terms.
        
        
        return: list, pauli sum equal to the input with like terms collected and sorted.
    
    '''
    simplified_terms = {}
    
    for pauli_string, coeff in terms:
        if pauli_string in simplified_terms:
            simplified_terms[pauli_string] += coeff
        else:
            simplified_terms[pauli_string] = coeff
    
    simplified_terms_list = [(pauli_string, coeff) for pauli_string, coeff in simplified_terms.items()]
    
    # will sort alphabetically
    return sorted(simplified_terms_list)



def test_jw(is_debug=True):
    # get the integrals

    geometry = "H 0 0 0; H 0 0 0.7414"
    basis = "sto-3g"
    charge = 0

    mol = pyscf.gto.Mole()
    # mol.unit = "bohr" # Default is angstrom
    mol.build(atom=geometry, basis=basis, charge=charge)

    # h, u = get_hs(mol, is_rhf=True)
    h, u = get_h2()


    # print(h.shape, u.shape)


    # ------------test my jordan_wigner---------------

    h = Hamiltonian(h, u)
    jw_h = jordan_wigner(h)

    if is_debug:
        print('My jordan wigner')
    for op, coeff in jw_h:
        if is_debug:
            print(f'{coeff:+.8f} * {op}')

    # get the transformation with qiskit

    
    '''From qiskit, for testing the jw'''
    driver = PySCFDriver(atom="H 0 0 0; H 0 0 0.7414", charge=0, spin=0, basis='sto3g')
    problem = driver.run()

    hamiltonian = problem.hamiltonian.second_q_op()
    
    # for label, coeff in sorted(hamiltonian.items()):
        # print(f"{coeff:+.8f} * '{label}'")


    mapper = JordanWignerMapper()
    qubit_op = mapper.map(hamiltonian)
    H2_pauli = []

    if is_debug:
        print('qiskit')
    for pauli, coeff in sorted(qubit_op.label_iter()):
        if is_debug:
            print(f"{coeff.real:+.8f} * {pauli}")
        H2_pauli.append((pauli, coeff))

    print(f"my_pauli_len={len(jw_h)}, qiskit_len={len(H2_pauli)}")

if __name__ == '__main__':
    import numpy as np

    from qiskit_nature.second_q.drivers import PySCFDriver
    from qiskit_nature.second_q.mappers import JordanWignerMapper
    import pyscf

    from pyscf_mel import get_hs
    from __qiskit_hamiltonian import get_h2
    

    test_jw(is_debug=True)
    

    