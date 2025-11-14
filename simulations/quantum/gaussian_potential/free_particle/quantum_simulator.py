from qiskit import QuantumCircuit as QC
from qiskit import transpile
from qiskit.circuit.library import iSwapGate
import numpy as np
from itertools import product
from qiskit_aer import QasmSimulator

class QuantumSimulator:
    def __init__(self, qc: QC, shots: int, psi: np.ndarray):
        self.qc = qc
        self.psi = psi
        self.shots = shots
        self.backend = QasmSimulator(method='statevector')

    def run(self):
        job = self.backend.run(self.qc, shots=self.shots)
        result = job.result()
        counts = result.get_counts()
        return counts
    
    def post_process(self, counts):
        """Post-process the measurement counts to extract non-zero amplitudes"""
        states = [''.join(bits) for bits in product('01', repeat=self.qc.num_qubits)]
        probabilities = np.array([counts.get(s, 0) / self.shots for s in states])
        states_arr = np.array(states)

        mask = probabilities > 0
        pvec_nonzero = probabilities[mask]
        states_nonzero = states_arr[mask]

        amps_nonzero = np.sqrt(pvec_nonzero)

        # --- pretty print results
        # print("Non-zero states:", states_nonzero)
        # print("Probabilities:", pvec_nonzero)
        # print("Amplitudes:", amps_nonzero)

        # --- optional: combine into a dictionary
        filtered = dict(zip(states_nonzero, pvec_nonzero))
        filtered_amps = dict(zip(states_nonzero, amps_nonzero))
        # print("Dictionary (probabilities):", filtered)
        # print("Dictionary (amplitudes):", filtered_amps)

        return states_nonzero, amps_nonzero
    
    def dispersion_factor(self, final_amps: np.ndarray):
        """
        Compare initial and final amplitude distributions to quantify dispersion.
        Returns a number in [0, 1]:
            0 = identical distributions (no spread)
            1 = fully different / maximally dispersed
        """
        # --- Normalize both wavefunctions
        p_init = np.abs(self.psi)**2
        p_init /= p_init.sum()
        
        p_final = np.abs(final_amps)**2
        p_final /= p_final.sum()
        
        # --- Compute cosine similarity of the two probability vectors
        numerator = np.sum(p_init * p_final)
        denom = np.sqrt(np.sum(p_init**2) * np.sum(p_final**2))
        similarity = numerator / denom if denom != 0 else 0.0
        
        # --- Dispersion factor (1 - similarity)
        dispersion = 1 - similarity
        return dispersion