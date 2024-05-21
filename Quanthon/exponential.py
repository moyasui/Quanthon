import numpy as np

def exponential_pauli(qc, pauli_str, coeff, method='inverted staircase'):
	"""The exponential of an pauli string, the circuit of qc is updated in place
	args:
		qc: Qubits object, the Qubits object the exponential is acting on,
		pauli_str: str, the Pauli string to exponentiate,
		coeff: float, the coefficient of the Pauli string,
		method: str, the method used to exponentiate the pauli string, defaults to 'inverted staircase'.
		
	"""
	if method == 'staircase':
		staircase(qc, pauli_str, coeff)
		
	elif method == 'inverted staircase':
		inverted_staircase(qc, pauli_str, coeff)

	elif method == 'fswap':
		pass
	else:
		raise ValueError('Invalid method, must be one of "staircase", "inverted staircase", or "fswap".')

def staircase(qc, pauli_str, coeff):
	# left 
	for i, p in enumerate(pauli_str):
		if p == 'X':
			qc.H(i)
		elif p == 'Y':
			qc.Rz(-np.pi/2, i)
			qc.H(i)

	for i in range(1, qc.n_qubit):
		if pauli_str[i] == 'I':
			qc.SWAP(i-1, i)
		else:
			qc.CNOT(i-1, i)

	# Rz
	qc.Rz(2 * coeff, qc.n_qubit - 1)

	# right
	for i in range(qc.n_qubit - 1, 0, -1):
		if pauli_str[i] == 'I':
			qc.SWAP(i-1, i)
		else:
			qc.CNOT(i-1, i)
	
	for i, p in enumerate(pauli_str):
		if p == 'X':
			qc.H(i)
		elif p == 'Y':
			qc.H(i)
			qc.Rz(np.pi/2, i)
	

def inverted_staircase(qc, pauli_str, coeff):

	# left 
	for i, p in enumerate(pauli_str):
		if p == 'Z':
			qc.H(i)
		elif p == 'Y':
			qc.H(i)
			qc.Rz(-np.pi/2, i)

	for i in range(1, qc.n_qubit):
		if pauli_str[i] == 'I':
			qc.SWAP(i, i-1)
		else:
			qc.CNOT(i, i-1)

	# Rz
	qc.Rx(2 * coeff, qc.n_qubit - 1)

	# right
	for i in range(qc.n_qubit - 1, 0, -1):
		if pauli_str[i] == 'I':
			qc.SWAP(i-1, i)
		else:
			qc.CNOT(i, i-1)
	
	for i, p in enumerate(pauli_str):
		if p == 'X':
			qc.H(i)
		elif p == 'Y':
			qc.H(i)
			qc.Rz(-np.pi/2, i)

if __name__ == '__main__':

	from Quanthon import Qubits
	from Quanthon.base import Gate
	from Quanthon.utils import get_pauli
	from scipy.linalg import expm

	n = 1
	qc = Qubits(n)
	pauli_str = 'Y'

	a = 0.5
	coeff = -1j * a
	exponential_pauli(qc, pauli_str, a, method='staircase')
	qc.run()
	print(qc)

	qc = Qubits(n)
	qc.reset_circuit()
	qc.circuit.append(Gate(f'exp({pauli_str})', expm(coeff * get_pauli(pauli_str)), n_qubits=qc.n_qubit))
	qc.run()
	print(qc)



