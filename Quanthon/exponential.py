import numpy as np

def exponential_pauli(qc, pauli_str, coeff, method='inverted staircase'):
	"""The exponential of an pauli string, the circuit of qc is updated in place. 
	Must be used for at least two qubits.
	args:
		qc: Qubits object, the Qubits object the exponential is acting on,
		pauli_str: str, the Pauli string to exponentiate,
		coeff: float, the coefficient of the Pauli string,
		method: str, the method used to exponentiate the pauli string, defaults to 'inverted staircase'.
		
	"""
	if qc.n_qubit < 2:
		raise ValueError('For the exponential of a single pauli operator, simply use one of the built-in gates Rx, Ry or Rz.')
	
	if method == 'staircase':
		staircase(qc, pauli_str, coeff)
		
	elif method == 'inverted staircase':
		inverted_staircase(qc, pauli_str, coeff)

	elif method == 'fswap':
		raise NotImplementedError("Not implemented.")
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
			qc.Rz(-np.pi/2, i)

	for i in range(1, qc.n_qubit):
		if pauli_str[i] == 'I':
			qc.SWAP(i-1, i)
		else:
			qc.CNOT(i, i-1)

	# Rx
	qc.Rx(2 * coeff, qc.n_qubit - 1)

	# right
	for i in range(qc.n_qubit - 1, 0, -1):
		if pauli_str[i] == 'I':
			qc.SWAP(i-1, i)
		else:
			qc.CNOT(i, i-1)
	
	for i, p in enumerate(pauli_str):
		if p == 'Z':
			qc.H(i)
		elif p == 'Y':
			qc.Rz(np.pi/2, i)





