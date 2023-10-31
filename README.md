# Quanthon

A minimal Python library for quantum computing, for physicists!

## Installation

Using PIP:
```sh
pip3 install Quanthon
```

## Get Started

### Importing the Module

```python
import Quanthon as qt
```

### Initializing a Single Qubit

Initialize a single qubit by creating an instance of the `Qubits` class.

```python
qubit = qt.Qubits(1)
```

### Applying Quantum Gates

`Quanthon` supports various quantum gates like the Hadamard (H), Pauli-X (X), Pauli-Y (Y), and Pauli-Z (Z) gates.

```python
# Apply a Hadamard gate on the first qubit
qubit.hadamard(0)

# Apply a Pauli-X gate on the first qubit
qubit.X(0)

# Apply a Pauli-Y gate on the first qubit
qubit.Y(0)

# Apply a Pauli-Z gate on the first qubit
qubit.Z(0)
```

### Performing Quantum Measurements

You can perform quantum measurements on your qubit system with a specific number of shots.

```python
result = qubit.measure(n_shots=10)
```

### Working with Multiple Qubits

Use the `Qubits` class for all states.

```python
from Quanthon import Qubits

# Initialize a 2-qubit system
two_qubits = Qubits_2()

# Initialize an n-qubit system
n_qubits = Qubits(3)
```

### CNOT and Swap Operations for Multiple Qubits

`Quanthon` allows you to perform CNOT and SWAP operations on multi-qubit systems.

```python
# Perform a CNOT operation between the first and second qubit
two_qubits.cnot(0, 1)

two_qubits.swap(0, 1)
```

