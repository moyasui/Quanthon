from itertools import product
import numpy as np


def get_all_pauli(n):

	'''Returns a list of all possible Pauli strings of length n (number of qubits)'''
	pauli = 'IXYZ'
	# all_paulis = list(product(pauli, pauli, repeat=int(n/2)))
	perms = [''.join(p) for p in product(pauli, repeat=n)]
	return perms


def _no_pauli_after_i(pauli_str):

	detected_I = False
	for op in pauli_str:
		if not detected_I:
			if op == 'I':
				detected_I = True
				continue
			else:
				continue
		if op != 'I':
			return False
	return True

def get_SWAPs(pauli_str, SWAP_list):

	if _no_pauli_after_i(pauli_str):
		return ''.join(pauli_str), SWAP_list

	found_I = False
	pauli_str = list(pauli_str)
	non_i_indx = None
	for i, op in enumerate(pauli_str):
		if op == 'I':
			found_I = True
			I_indx = i
			continue
		if found_I:
			if op != 'I':
				non_i_indx = i
			pauli_str[I_indx] = op 
			pauli_str[non_i_indx]= 'I'
			SWAP_list.append((I_indx, non_i_indx))
			found_I = False

	return get_SWAPs(pauli_str, SWAP_list)

def get_CNOTs(pauli_str): 

	if not _no_pauli_after_i:
		raise ValueError('Pauli string must not have any Pauli operator after I.')
	CNOT_pairs = []

	for i in range(len(pauli_str)-1):
		if pauli_str[i+1] == 'I':
			break
		CNOT_pairs.append((i+1, i))
		
	return CNOT_pairs

def change_basis(pauli_str):

	for i in pauli_str:
		if i not in 'IXYZ':
			raise ValueError(f'Invalid Pauli operator: {i}')

	SWAPed_pauli, SWAPs = get_SWAPs(pauli_str, [])
	CNOT_pairs = get_CNOTs(SWAPed_pauli)

	change_basis_op = []

	for pair in SWAPs:
		change_basis_op.append(('SWAP', pair))
	
	# print(SWAPed_pauli)
	for i, op in enumerate(SWAPed_pauli):
		if op == 'X':
			change_basis_op.append(('H', i))
		elif op == 'Y':
			change_basis_op.append(('Sdag', i))
			change_basis_op.append(('H', i))
		else:
			change_basis_op.append(('I', i))
 
	
	for pair in reversed(CNOT_pairs):
		change_basis_op.append(('CNOT', pair))
	
	return change_basis_op

def rotate_basis(qc, ops):
	'''rotate to the measurement basis'''
	
	# print("bc", sum(qc.state**2))

	for op_type, idx in ops:
		if op_type == 'CNOT':
			qc.CNOT(idx[0], idx[1])
		elif op_type == 'SWAP':
			qc.SWAP(idx[0], idx[1])
		elif op_type == 'X':
			qc.X(idx)
		elif op_type == 'Y':
			qc.Y(idx)
		elif op_type == 'Z':
			qc.Z(idx)
		elif op_type == 'H':
			qc.H(idx)
		elif op_type == 'Sdag':
			qc.Sdag(idx)
		elif op_type == 'I':
			continue
		else:
			raise ValueError(f'Invalid operator: {op_type}')
	
	qc.run()
	# print([(gate.name, gate.n_qubits) for gate in qc.circuit])
	# print('ops', ops)
	# print([gate.name for gate in qc.circuit])
	# for gate in qc.circuit:
	#     qc.state = gate.act(qc.state)
		# print(gate.name)
		# print(qc.state)
		# print(qc.state**2)
		# print(sum(qc.state**2))
	
	
	
def cal_expectation(qc, pauli_ops, n_shots=10000):
	
	'''
	return:
		the expectation value
	args:
		qc: the circuit to be measured in
	'''

	expectation = 0
	for pauli_str, coeff in pauli_ops:

		
		if set(pauli_str) == {'I'}:
			expectation += coeff
			continue
		# print(pauli_str, coeff)
		state_eigval = np.ones(len(qc.state))
		state_eigval[int(0.5*len(qc.state)):] *= -1 # first half of state has eigenvalue 1 and the second half -1
		# print(state_eigval)
		qc_copy = qc.copy()
		cb_ops = change_basis(pauli_str)
		# print(cb_ops)

		rotate_basis(qc_copy, cb_ops)

		# print(qc_copy)
		counts = qc_copy.measure(n_shots)[:, 0]
		# print("this count", counts)
		# print(counts)
		expectation += coeff * np.sum(counts * state_eigval) / n_shots
	
	# print(cb_ops)
	return expectation


if __name__ == '__main__':
	# Test Pauli operators
	
	# print(get_all_pauli(2))
	# pauli_str = 'YY'
	# for op in pauli_str:
	#     rotation_ops = get_rotations(op)
	#     print(rotation_ops)

	a = 'X'
	b = 'Y'
	c = 'IZ'
	d = 'ZIXX'
	# test no_op_after_i
	
	# print(_no_pauli_after_i(a), _no_pauli_after_i(b), _no_pauli_after_i(c))

	# print(_no_pauli_after_i(d))

	# PASSED

	# test get_SWAPs
	
	# for i in [a, b, c, d]:
	#     result, SWAPs = get_SWAPs(i, [])
	#     result = "".join(result)
	#     print(result, SWAPs)
	
	# PASSED

	# test get_CNOTs

	# for paulis in [a, b, c, d]:
	#     result, SWAPs = get_SWAPs(paulis, [])
	#     result = "".join(result) 
	#     print(result)
	#     CNOT_pairs = get_CNOTs(result)
	#     print(f"CNOTs: {paulis}", CNOT_pairs)

	# PASSED

	# test change_basis

	# for i in [a, b, c, d]:
	#     big_op = change_basis(i)
	#     print(f"big_op for {i}", big_op)

	# PASSED

	# test expectation

	# a = 'X'
	# b = 'IYYI'
	# c = 'ZIII'
	# d = 'XIXI'

	# b = 'XX'
	# c = 'YY'
	# d = 'IZ'
	# e = 'XY'

	# testing drawing, seems ok

	# qc = Qubits(2)
	# qc.H(0)
	# qc.CNOT(0,1)
	# qc.run()
	# # print(qc)

	# hamiltonian_str = [b, c, e]
	# coeffs = {b:0.5 , c:0.5 , d:0.5, e:3}
	# ham = [(op, coeffs[op]) for op in hamiltonian_str]
	
	# # print(ham)

	# print(f"hamiltonian: {hamiltonian_str}")


	# # print(ham_mat)
	# # for term in hamiltonian:
		
	# #     qc_copy = qc.copy()
	# #     print(qc_copy)
	# #     ops = change_basis(term)
	# #     rotate_basis(qc_copy, ops)
	# #     print(ops)
	# #     print(qc_copy.gate_history)
	# #     qc_copy.draw()

	# # testing expectation
	
	# for term in hamiltonian_str:
		
	#     qc_copy = qc.copy()
	#     # print(qc_copy)
	#     ops = change_basis(term)
	#     rotate_basis(qc_copy, ops)
	#     # print(ops)
	
	# estimated_energy = expectation(qc_copy, ham, n_shots=100000)
	# print(estimated_energy)

	# from Quanthon import pauli_sum
	# # ham_mat = pauli_sum(ham)
	# # real_energy = qc.state.conj() @ ham_mat @ qc.state
	# # print(real_energy.real)

	# # Testing expectation 2

	# qc = Qubits(4)
	# qc.H(0)
	# qc.CNOT(0,1)
	# qc.CNOT(1,2)
	# qc.CNOT(2,3)
	# qc.run()
	# Ham = [
	#     # ('ZIZI', (1)),
	#     ('ZZXX', (0.5)), 
	#     # ('YYII', (1)),
	#     # ('IXIX', (0.5))
	#     ]

	# H_mat = pauli_sum(Ham)
	# ee = expectation(qc, Ham, n_shots=100000)
	# re = qc.state.conj().T @ H_mat @ qc.state

	# print("exact: ", re, "expectation: ", ee)
 


