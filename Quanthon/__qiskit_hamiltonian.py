from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper

import numpy as np

def get_h2():

    driver = PySCFDriver(atom="H 0 0 0; H 0 0 0.7414", charge=0, spin=0, basis='sto-3g')
    problem = driver.run()

    hamiltonian = problem.hamiltonian.second_q_op()

    h = np.zeros((4,4))
    u = np.zeros(((4,4,4,4)))
    for label, coeff in sorted(hamiltonian.items()):
        label = remove_non_num(label)
        if len(label) == 2:
            h[int(label[0]), int(label[1])] = coeff
        elif len(label) == 4:
            u[int(label[0]), int(label[1]), int(label[2]), int(label[3])] = coeff
        else:
            raise ValueError("Label not recognized")
        # print(f"{coeff:+.8f} * '{label}'")

    return h, u

def remove_non_num(string):

    for char in string:
        if char not in '0123456789':
            string = string.replace(char, '')
    return string


if __name__ == "__main__":
    h, u  = get_h2()

    u_list = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    # if u[i,j,k,l] != 0:
                    #     print(f'{i}{j}{k}{l}, {u[i,j,k,l]}')
                    if len({i,j,k,l}) == 1:
                        continue
                    
                    if (i == k and j == l) or (i == l and j == k):
                        print(i, j, k, l)
                        u_list.append(u[i,j,k,l])

    u_list = 0.25 * np.array(u_list)
    print(sum(u_list))
    # get this -0.81261796

    h_coeff = h.diagonal()
    h_coeff = 0.5 * h_coeff
    
    print(sum(h_coeff))

    print(sum(u_list) + sum(h_coeff))