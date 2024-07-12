'''Doc:
	Ansatz for VQE calculations. Can grow if needs be. Parametrise the qubits object.
	parameters are taken as arguments to the circuit at every iteration'''

from .mappers import jordan_wigner
from .physics import Hamiltonian

import numpy as np
from scipy.linalg import expm


from .base import Qubits, Gate
from .utils import pauli_sum
from .exponential import exponential_pauli

# DEBUGGER

def qprint(msg):
	
	print("QubitAdaptAnsatz")
	print(msg)


class Ansatz:

	'''Should work for any type of ansatz that are not evolving, do not use on its own.'''
	def __init__(self, n_qubits, n_params, reps=1) -> None:

		self.n_qubits = n_qubits
		self.reps = reps
		self.n_params = n_params

	def run(self):
		self.qubits.run()
		
class HardwareEfficientAnsatz(Ansatz):

	def __init__(self, n_qubits, reps=1) -> None:
		super().__init__(n_qubits, 2 * n_qubits * reps, reps)


	def create_circuit(self, init_params, init_state=None):
		

		self.qubits = Qubits(self.n_qubits)
		if init_state is not None:
			self.qubits.set_state(init_state)
		if len(init_params) != self.n_params:
			raise ValueError(f'Initial parameters do not match the number of parameters required for the ansatz: {len(init_params)} != {self.n_params}')
		
		reshaped_params = init_params.reshape(self.reps, 2*self.n_qubits)
		
		for r in range(self.reps):
			for i in range(self.n_qubits):
				# check if the parameters are real
				# if reshaped_params[r, i] != reshaped_params[r, i].real:
				# 	raise ValueError(f'Parameter {[r,i]} is not real.')
				self.qubits.Rx(reshaped_params[r, i], i)
				self.qubits.Ry(reshaped_params[r, i + self.n_qubits], i)
		
		for r in range(self.reps):
			for i in range(self.n_qubits):
				if i != self.n_qubits - 1:
					self.qubits.CNOT(i, i+1)

		# self.qubits.draw()
		
class RyAnsatz(Ansatz):
		
	def __init__(self, n_qubits, reps=1) -> None:
		super().__init__(n_qubits, n_qubits*reps, reps)
	
	def create_circuit(self, init_params, init_state=None):

		if init_state is not None:
			self.qubits.set_state(init_state)
		
		self.qubits = Qubits(self.n_qubits)
		if len(init_params) != self.n_params:
			raise ValueError(f'Initial parameters do not match the number of parameters required for the ansatz: {len(init_params)} != {self.n_params}')
		
		reshaped_params = init_params.reshape(self.reps, self.n_qubits)
		
		for r in range(self.reps):
			for i in range(self.n_qubits):
				self.qubits.Ry(reshaped_params[r, i], i)

		for r in range(self.reps):
			for i in range(self.n_qubits):
				if i != self.n_qubits - 1:
					self.qubits.CNOT(i, i+1)
			# self.qubits.draw()
			
class UCCSDAnsatz(Ansatz):

	def __init__(self, n_qubits, n_params, n_occupied, n_virtual, reps=1,) -> None:
		super().__init__(n_qubits, n_params, reps)
		self.n_occupied = n_occupied
		self.n_virtual = n_virtual
	
	def get_singles(self):

		coeffs = np.ones((self.n_occupied, self.n_virtual))
		tb = np.zeros((self.n_occupied, self.n_virtual, self.n_occupied, self.n_virtual))
		h = Hamiltonian(coeffs, tb)
		print(h)
		singles = jordan_wigner(h)
		print(singles)

	
		return singles

	def get_doubles(self):

		coeffs = np.ones((self.n_occupied, self.n_virtual))
		tb = np.zeros((self.n_occupied, self.n_virtual, self.n_occupied, self.n_virtual))
		h = Hamiltonian(coeffs, tb)
		print(h)
		singles = jordan_wigner(h)
		print(singles)
		
		return doubles
		
	def create_circuit(self, init_params, init_state=None):
		if init_state is not None:
			self.qubits.set_state(init_state)
		
		self.qubits = Qubits(self.n_qubits)
		if len(init_params) != self.n_params:
			raise ValueError(f'Initial parameters do not match the number of parameters required for the ansatz: {len(init_params)} != {self.n_params}')
		
		reshaped_params = init_params.reshape(self.reps, self.n_qubits)
		
		# singles






class QubitAdaptAnsatz:

	def __init__(self, n_qubits, pool='V', init_state=None) -> None:

		'''Ansatz for the Adapt-VQE calculation.
		args:
			n_qubits: int, number of qubits in the system
			pool: str, the pool of operators to choose from, default is 'V' pool which can't be used fewer than 3 qubits!!.
			init_state: np.array, the initial state of the system, default is None
		'''

		self.qubits = Qubits(n_qubits)

		if init_state is None:
			# initial state not sepcified, initialise randomly
			# init_state = np.random.rand(2**n_qubits) + 1j * np.random.rand(2**n_qubits)
			init_state = np.ones(2**n_qubits)
			init_state = init_state / np.linalg.norm(init_state)
		
		self.init_state = init_state
		self.qubits.set_state(init_state)

		if pool == 'V':
			self.pool = self.create_complete_V_pool(n_qubits)
		
		elif pool == 'Vx':
			self.pool = self.create_complete_Vx_pool(n_qubits)
			
		elif pool == 'G':
			self.pool = self.create_complete_G_pool(n_qubits)
		
		else:
			self.pool = pool # custom pool, doesn't have to be complete
			
		

	def __repr__(self) -> str:
		
		return f"QubitAdaptAnsatz of {self.qubits} qubits and pool: {self.pool}."
	
	def create_complete_V_pool(self, n):
		
		if n == 1:
			return ['iY', 'Z'] 
		
		if n == 2:
			return ['iYZ', 'iIY']

		prev_set = self.create_complete_V_pool(n - 1)
		new_pool = [] # set()

		for prev_op in prev_set:
			new_op = prev_op + 'Z'
			new_pool.append(new_op)


		new_pool.append('i' + (n-1) * 'I' + 'Y')
		new_pool.append('i' + (n-2) * 'I' + 'YI')


		return new_pool
		
	def create_complete_Vx_pool(self, n):
		
		if n == 1:
			return ['iY', 'X'] 
		
		if n == 2:
			return ['iYX', 'iIY']

		prev_set = self.create_complete_Vx_pool(n - 1)
		new_pool = [] # set()

		for prev_op in prev_set:
			new_op = prev_op + 'X'
			new_pool.append(new_op)


		new_pool.append('i' + (n-1) * 'I' + 'Y')
		new_pool.append('i' + (n-2) * 'I' + 'YI')


		return new_pool
	
	def create_complete_G_pool(self, n):
		
		if n == 1:
			return ['iY', 'Z']
		pool = self._get_ys(n)
		# print("get ys", pool)
		pool.extend(self._get_yzs(n))
		# print("get yzs", pool)
		return pool

	
	def _get_ys(self,n):
		'''Get a list of all the operators with a single Y in them'''
		op = n*'I'
		ops = []
		for i in range(1, n):
			y_op = list(op)
			y_op[i] = 'Y'
			ops.append('i' + ''.join(y_op))
		
		return ops

	def _get_yzs(self, n):
		op = n*'I'
		ops = []
		for i in range(n-1):
			yz_op = list(op)
			yz_op[i] = 'Y'
			yz_op[i+1] = 'Z'
			ops.append('i' + ''.join(yz_op))
		return ops

	def append_op(self, op, decompose_exp, decompose_method, print_log=False):

		'''
		op: string, representing one of the operators in the pool
		
		'''
		# no longer need to append all the gates since they are now saved

		if print_log:
			qprint(f"append_op op: {op}")

		if decompose_exp:
			exponential_pauli(self.qubits, op.strip('i'), coeff=None, method=decompose_method)
			if print_log:
				qprint(f"len circuit: {len(self.qubits.circuit)}")
		else:	
			op_mat = pauli_sum([(op.strip('i'), 1j)])
			def parametrised_mat(param, op_mat):
				return expm(-0.5*param * op_mat)
			
			adapt_gate = Gate(f'exp({op})', matrix=lambda param: parametrised_mat(param, op_mat), n_qubits=self.qubits.n_qubit, is_parametrised=True)
			self.qubits.circuit.append(adapt_gate)
			
	
	def run(self, params):
		param_indx = 0
		for gate in self.qubits.circuit:
			# print(gate.matrix(params[i]))
			if gate.is_parametrised:

				self.qubits.state = gate.act(self.qubits.state, params[param_indx])
				param_indx += 1
			else:
				self.qubits.state = gate.act(self.qubits.state)
		
		# print(self.qubits.circuit)
	
	def run_without_update(self, params):
		'''Run the circuit without updating the state.'''
		state = self.qubits.state
		param_indx = 0

		for gate in self.qubits.circuit:
			if gate.is_parametrised:
				state = gate.act(state, params[param_indx])
				param_indx += 1
			else:
				state = gate.act(state)
			
			if not np.isclose(np.linalg.norm(state), 1, atol=1e-4):
				raise ValueError(f"Probability does not sum to unity, prob: {self.qubits.prob}, state: {self.qubits.state}.")
		
		return state



	


	




