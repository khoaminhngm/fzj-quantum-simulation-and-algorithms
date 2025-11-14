import numpy as np
from qiskit import QuantumCircuit as QC


class TimeEvolution:
    def __init__(self, qc: QC, psi: np.ndarray, b: float, a: float):
        self.psi = psi
        self.a = a
        self.b = b
        self.num_qubits = int(np.log2(len(psi)))
        self.qc = qc

    def _time_transform(self):
        theta = 2 * np.arctan2(np.abs(self.a), np.abs(self.b))
        # theta = np.pi/2  # Example fixed angle for demonstration with iSwap
        n = range(self.qc.num_qubits - 1)
        N = self.qc.num_qubits
        for i in n:
            self.qc.rxx(-theta, self.qc.qubits[i], self.qc.qubits[i+1])
            self.qc.ryy(-theta, self.qc.qubits[i], self.qc.qubits[i+1])
        
        self.qc.rxx(-theta, self.qc.qubits[N-1], 0)
        self.qc.ryy(-theta, self.qc.qubits[N-1], 0)
        return self.qc

    def apply_time_evolution(self, time_steps: int) -> QC:
        for i in range(time_steps):
            self.qc = self._time_transform()
        return self.qc