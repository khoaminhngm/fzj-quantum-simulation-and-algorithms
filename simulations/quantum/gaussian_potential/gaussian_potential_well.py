import numpy as np
from qiskit import QuantumCircuit

class GaussianPotential:
    def __init__(self, V0, sigma, x0, wave_function: np.ndarray, qc: QuantumCircuit):
        self.V0 = V0
        self.sigma = sigma
        self.x0 = x0
        self.wave_function = wave_function
        self.qc = qc
        print("Gaussian Potential Well Initialized")
    
    def _potential(self, x):
        return -self.V0 * np.exp(-((x - self.x0) ** 2) / (2 * self.sigma ** 2))
    
    def create_potential_matrix(self):
        potential_values = [self._potential(x) for x in range(len(self.wave_function))]
        U = np.diag(potential_values)
        print("Gaussian Potential Matrix Created")
        return U
    
    def apply_potential(self, dt):
        for i in range(len(self.wave_function)):
            U_i = self._potential(i)
            self.qc.rz(-dt * U_i, i)  # minus from derivation (RZ(-α))
        return self.qc


if __name__ == "__main__":
    V0 = 10.0
    sigma = 3
    psi = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.7, 0.9, 0.7, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1], dtype=complex)
    psi /= np.linalg.norm(psi)
    x0 = len(psi) // 2  

    qc = QuantumCircuit(len(psi))
    gaussian_potential = GaussianPotential(V0, sigma, x0, psi, qc)

    U = gaussian_potential.create_potential_matrix()
    print("\nPotential Matrix:\n", U)

    qc = gaussian_potential.apply_potential(dt=0.1)
    print("\nCircuit:\n")
    print(qc)
