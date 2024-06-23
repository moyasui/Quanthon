import numpy as np
from Quanthon.base import Gate

def no_i_after_pauli(pauli_str):

	detected_p = False
	for op in pauli_str:
		if op != 'I':
			detected_p = True
			continue

		if detected_p and op == 'I':
			return False
	return True

def get_SWAPs(pauli_str, SWAP_list):

	if no_i_after_pauli(pauli_str):
		return ''.join(pauli_str), SWAP_list

	found_p = False
	pauli_str = list(pauli_str)
	for i, op in enumerate(pauli_str):
		if op != 'I':
			found_p = True
			continue
		if op == 'I' and found_p:
			SWAP_list.append((i, i-1))
			pauli_str[i] = pauli_str[i-1]
			pauli_str[i-1] = 'I'
			break

	return get_SWAPs(pauli_str, SWAP_list)

def get_CNOTs(pauli_str, direction='normal'): 

	if not no_i_after_pauli(pauli_str):
		raise ValueError('Pauli string must not have any Pauli operator after I.')
	CNOT_pairs = []

	if direction == 'normal':
		for i in range(len(pauli_str)-1):
			if pauli_str[i] == 'I':
				continue
			CNOT_pairs.append((i, i+1))
	elif direction == 'reverse':
		for i in range(len(pauli_str)-1, 0, -1):
			if pauli_str[i-1] == 'I':
				break
			CNOT_pairs.append((i, i-1))
			
	return CNOT_pairs

def exponential_pauli(qc, pauli_str, coeff=None, method='inverted staircase'):
	"""The exponential of an pauli string, the circuit of qc is updated in place. 
	Must be used for at least two qubits.
	args:
		qc: Qubits object, the Qubits object the exponential is acting on,
		pauli_str: str, the Pauli string to exponentiate,
		coeff: float, the coefficient of the Pauli string, a in e^{-iaP}
		method: str, the method used to exponentiate the pauli string, defaults to 'inverted staircase'.
		
	"""
	if qc.n_qubit < 2:
		if pauli_str == 'I':
			if coeff is None:
				gate_mat = lambda coeff: (np.cos(coeff) - 1j * np.sin(coeff)) * np.eye(2)
				qc.circuit.append(Gate(f'exp(i{pauli_str})', gate_mat, n_qubits=qc.n_qubit, is_parametrised=True))
			else:
				qc.circuit.append(Gate(f'exp(i{pauli_str})', np.exp(-1j * coeff) * np.eye(2), n_qubits=qc.n_qubit))
		elif pauli_str == 'X':
			if coeff is None:
				qc.Rx(None, 0)
			else:
				qc.Rx(2 * coeff, 0)
		elif pauli_str == 'Y':
			if coeff is None:
				qc.Ry(None, 0)
			else:
				qc.Ry(2 * coeff, 0)
		elif pauli_str == 'Z':
			if coeff is None:
				qc.Rz(None, 0)
			else:
				qc.Rz(2 * coeff, 0)

		return

	if method == 'staircase':
		staircase(qc, pauli_str, coeff)
		
	elif method == 'inverted staircase':
		inverted_staircase(qc, pauli_str, coeff)
	

	else:
		raise ValueError('Invalid method, must be one of "staircase", or "inverted staircase".')


def staircase(qc, pauli_str, coeff):

	# left 

	new_ps, swaps = get_SWAPs(pauli_str, [])
	CNOT_pairs = get_CNOTs(new_ps, direction='normal')
	
	for pair in swaps:
		qc.SWAP(pair[0], pair[1])

	for i, p in enumerate(new_ps):
		if p == 'X':
			qc.H(i)
		elif p == 'Y':
			qc.Rz(-0.5 * np.pi, i)
			qc.H(i)

	for cnot in CNOT_pairs:
		qc.CNOT(cnot[0], cnot[1])

	# Rz
	if coeff is None:
		qc.Rz(None, qc.n_qubit - 1)
	else:
		qc.Rz(2 * coeff, qc.n_qubit - 1)

	# right
	for cnot in reversed(CNOT_pairs):
		qc.CNOT(cnot[0], cnot[1])

	for i in range(qc.n_qubit - 1, -1, -1):
		if new_ps[i] == 'X':
			qc.H(i)
		elif new_ps[i] == 'Y':
			qc.H(i)
			qc.Rz(0.5 * np.pi, i)
	
	for pair in reversed(swaps):
		qc.SWAP(pair[0], pair[1])
	

def inverted_staircase(qc, pauli_str, coeff):
	
	# left 
	new_ps, swaps = get_SWAPs(pauli_str, [])
	CNOT_pairs = get_CNOTs(new_ps, direction='reverse')
	
	for pair in swaps:
		qc.SWAP(pair[0], pair[1])

	for i, p in enumerate(new_ps):
		if p == 'Z':
			qc.H(i)
		elif p == 'Y':
			qc.Rz(-0.5 * np.pi, i)

	for cnot in reversed(CNOT_pairs):
		qc.CNOT(cnot[0], cnot[1])

	# Rx
	if coeff is None:
		qc.Rx(None, qc.n_qubit - 1)
	else:
		qc.Rx(2 * coeff, qc.n_qubit - 1)

	# right
	for cnot in CNOT_pairs:
		qc.CNOT(cnot[0], cnot[1])

	for i in range(qc.n_qubit - 1, -1, -1):
		if new_ps[i] == 'Z':
			qc.H(i)
		elif new_ps[i] == 'Y':
			qc.Rz(0.5 * np.pi, i)
	
	for pair in reversed(swaps):
		qc.SWAP(pair[0], pair[1])

	


if __name__ == '__main__':

	from Quanthon import Qubits
	from Quanthon.base import Gate
	from Quanthon.utils import get_pauli
	from scipy.linalg import expm

	# paulis = 'IXYZIYIYZ'
	# new_ps, swaps = get_SWAPs(paulis, [])
	# print(get_CNOTs(new_ps, direction='reverse'))
	
	pauli_str = 'iXYZI'
	n = len(pauli_str) - 1
	pauli_str = pauli_str.strip('i')
	print(pauli_str)

	a = np.pi/3
	coeff = -1j * a


	state = np.ones(2**n)
	state /= np.linalg.norm(state)
	
	qc = Qubits(n)
	qc.set_state(state)
	exponential_pauli(qc, pauli_str, a, method='staircase')
	qc.run()
	print("staircase", qc)
	for gate in qc.circuit:
		print(gate.name)
		# print(qc)
		# qc.state = gate.act(qc.state)
	s_state = qc.state

	# with inverted staircase algorithm
	qc = Qubits(n)
	qc.set_state(state)
	exponential_pauli(qc, pauli_str, a, method='inverted staircase')
	# for gate in qc.circuit:
	# 	print(gate.name)
		# print(qc)
	qc.run()
	print("inverted", qc)
	is_state = qc.state

	# with scipy.linalg.expm
	qc = Qubits(n)
	qc.set_state(state)
	# print(qc)
	qc.circuit.append(Gate(f'exp(i{pauli_str})', expm(coeff * get_pauli(pauli_str)), n_qubits=qc.n_qubit))
	# for gate in qc.circuit:
	# 	print(gate)
		# qc.state = gate.act(qc.state)
	qc.run()
	print("scipy", qc)
	sp_state = qc.state

	# print(f"{sp_state-is_state}")
	# print(f"{sp_state-s_state}")

	print("invert", np.allclose(sp_state, is_state))
	print("stair", np.allclose(sp_state, s_state))


	# # parametrised exponential
	# qc = Qubits(n)
	# exponential_pauli(qc, pauli_str, None, method='inverted staircase')
	# # exponential_pauli(qc, pauli_str.strip('i'), coeff=None, method='inverted staircase')
	# qc.run([2*a, 2*a])
	# print("parametrised", qc)
