import numpy as np
from qiskit import QuantumCircuit as QC


class TimeEvolution:
    def __init__(self, qc: QC, psi: np.ndarray, b: float, a: float, time_steps: int):
        self.psi = psi
        self.U = np.matrix([[1, 0, 0, 0],
                             [0, b, a, 0],
                             [0, a, b, 1],
                             [0, 0, 0, 1]])
        self.time_steps = time_steps
        self.num_qubits = int(np.log2(len(psi)))
        self.qc = qc

    def time_transform(self):
        theta = np.arccos(self.U[1,1].real)
        # theta = np.pi/2  # Example fixed angle for demonstration with iSwap
        n = range(self.qc.num_qubits - 1)
        for i in n:
            self.qc.rxx(-theta, self.qc.qubits[i], self.qc.qubits[i+1])
            self.qc.ryy(-theta, self.qc.qubits[i], self.qc.qubits[i+1])
        return self.qc

    def apply_time_evolution(self) -> QC:
        for i in range(self.time_steps):
            self.qc = self.time_transform()
        return self.qc