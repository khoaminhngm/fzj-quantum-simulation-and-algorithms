from qiskit import QuantumCircuit as QC
import numpy as np

class State_Preparator:
    def __init__(self, psi: np.array, qc: QC):
        self.psi = psi / np.linalg.norm(psi) # normalize
        self.qc = qc

    def flip_qubits(self, qubits: list[int]) -> QC:
        for q in qubits:
            # print(f"Flipping q{q}")
            self.qc.x(q)
        return self.qc
    
    def prepare(self):
        controls = [0]
    
        # Init first qubit
        theta1 = 2*np.arcsin(np.abs(self.psi[0])).real
        self.qc.ry(theta1, 0)

        phi1 = np.angle(self.psi[0]) 
        self.qc.rz(phi1, 0)

        # # Prepare second qubit
        theta2 = 2*np.arcsin(np.abs(self.psi[1])/np.sqrt(1-np.abs(self.psi[0])**2)).real
        self.qc.cry(theta2, 0, 1, ctrl_state=0)

        phi2 = np.angle(self.psi[1])
        self.qc.crz(phi2, 0, 1, ctrl_state=0)

        # Subsequent qubits
        for i in range(1, self.qc.num_qubits - 1): # i+1 -> target
            controls.append(i)
            print("Controls:", controls)

            theta = 2*np.arcsin(
                self.psi[i+1]/np.sqrt(1-sum([np.abs(self.psi[j])**2 for j in range(i+1)]))
                ).real
            
            phi = np.angle(self.psi[i+1])
            
            # print("Theta-loop:", theta)

            # theta2 = 2*np.arcsin(self.psi[1]/np.sqrt(1-self.psi[0]**2)).real
            # print("Theta 2-manual:", theta2)
            # theta3 = 2*np.arcsin(self.psi[2]/np.sqrt(1-self.psi[0]**2-self.psi[1]**2)).real
            # print("Theta 3-manual:", theta3)

            self.flip_qubits(controls)
            self.qc.mcry(theta, controls, self.qc.qubits[i+1])
            self.qc.mcrz(phi, controls, self.qc.qubits[i+1])
            print(f"State prep for qubit {i+1} done.")
            self.flip_qubits(controls)

        return self.qc
    

if __name__ == "__main__":
    from measure import Measurement
    from quantum_simulator import QuantumSimulator

    qc = QC(3, 3)
    sp = State_Preparator(np.array([0.5, 0.7, 0.5]), qc)
    qc = sp.prepare()

    meas = Measurement(qc)
    qc = meas.measure_all()

    qs = QuantumSimulator(qc, shots=np.pow(2, 20), psi=np.array([0.5, 0.7j, 0.5]))
    counts = qs.run()
    states_nonzero, amps_nonzero = qs.post_process(counts)
    print("Non-zero states:", states_nonzero)
    print("Amplitudes:", amps_nonzero)