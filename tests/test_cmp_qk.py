from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.algorithms import GroundStateEigensolver

from qiskit.quantum_info import SparsePauliOp
from qiskit.primitives import Estimator
from qiskit.circuit.library import TwoLocal

from qiskit_algorithms import VQE
from qiskit_algorithms import AdaptVQE

from qiskit_algorithms.optimizers import COBYLA, L_BFGS_B, SLSQP
from qiskit_algorithms.utils import algorithm_globals

import warnings

import numpy as np

def qiskit_hardware(op):
        
        warnings.filterwarnings("ignore")

        estimator = Estimator()
        qubit_op = SparsePauliOp.from_list(op)

        # we will iterate over these different optimizers
        optimizers = [COBYLA(maxiter=80), L_BFGS_B(maxiter=60), SLSQP(maxiter=60)]
        converge_counts = np.empty([len(optimizers)], dtype=object)
        converge_vals = np.empty([len(optimizers)], dtype=object)

        for i, optimizer in enumerate(optimizers):
            print("\rOptimizer: {}        ".format(type(optimizer).__name__), end="")
            algorithm_globals.random_seed = 50
            ansatz = TwoLocal(rotation_blocks="ry", entanglement_blocks="cz")

            counts = []
            values = []

            def store_intermediate_result(eval_count, parameters, mean, std):
                counts.append(eval_count)
                values.append(mean)

            vqe = VQE(estimator, ansatz, optimizer, callback=store_intermediate_result)
            result = vqe.compute_minimum_eigenvalue(operator=qubit_op)
            converge_counts[i] = np.asarray(counts)
            converge_vals[i] = np.asarray(values)
        
        return result

def _qiskit_adapt(self):
    
    pass
    warnings.filterwarnings("ignore")

    estimator = Estimator()
    qubit_op = SparsePauliOp.from_list(
    [
        ("II", -1.052373245772859),
        ("IZ", 0.39793742484318045),
        ("ZI", -0.39793742484318045),
        ("ZZ", -0.01128010425623538),
        ("XX", 0.18093119978423156),
    ]
    )

    # we will iterate over these different optimizers
    optimizer = SLSQP(maxiter=60)

    print("\rOptimizer: {}        ".format(type(optimizer).__name__), end="")
    algorithm_globals.random_seed = 50
    ansatz = TwoLocal(rotation_blocks="ry", entanglement_blocks="cz")

    vqe = VQE(estimator, ansatz, optimizer)
    adapt_vqe = AdaptVQE(vqe)
    eigenvalue, _ = adapt_vqe.compute_minimum_eigenvalue(qubit_op)
    print(eigenvalue)

    return eigenvalue

def matrix_to_fermionic_op(one_body_matrix, two_body_matrix):
    # Get the size of the one-body matrix
    n = one_body_matrix.shape[0]

    # Initialize an empty FermionicOp dictionary
    op_dict = {}

    # Handle one-body terms
    for i in range(n):
        for j in range(n):
            if one_body_matrix[i, j] != 0:
                key = f"+_{i} -_{j}"
                op_dict[key] = one_body_matrix[i, j]

    # Handle two-body terms
    # Assuming two_body_matrix is a 4-dimensional array
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    if two_body_matrix[i, j, k, l] != 0:
                        # Create the corresponding fermionic operator for the two-body term
                        key = f"+_{i} +_{j} -_{k} -_{l}"
                        op_dict[key] = two_body_matrix[i, j, k, l]

    # Create and return the FermionicOp
    return FermionicOp(op_dict)

def make_random_hamiltonian(n=4, seed=42):

    # one_body = np.zeros((n,n))
    rng = np.random.default_rng(seed)
    one_body = rng.random((n,n)) + 1j*rng.random((n,n))
    one_body = one_body + one_body.conj().T # make it hermitian

    two_body = rng.random((n,n,n,n)) + 1j*rng.random((n,n,n,n))
    two_body = two_body + two_body.conj().T

    op = matrix_to_fermionic_op(one_body,two_body)

    h = Hamiltonian(np.flip(one_body), np.flip(two_body))

    # h: Hamiltonian, op: FermionicOp for qiskit stuff
    return h, op