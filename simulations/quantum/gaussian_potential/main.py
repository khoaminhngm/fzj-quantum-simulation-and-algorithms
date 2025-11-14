from free_particle import State_Preparator, TimeEvolution, Measurement, QuantumSimulator
from gaussian_potential_well import GaussianPotential
from plotter import WavePlotter
import numpy as np
from qiskit import QuantumCircuit
import matplotlib.pyplot as plt

def main():
    psi = np.array([0.5, 0.7, 0.5], dtype=complex) # complex wavefunction
    psi /= np.linalg.norm(psi)  # Normalize the wavefunction
    print(psi)
    b = 1/np.sqrt(2) # probability that the particle stays in the same position
    a = 1j/np.sqrt(2) # probability that the particle moves to the left/right
    time_steps = 2
    shots = np.pow(2, 13)

    # Initialize Classes
    # qc = QuantumCircuit(len(psi), len(psi))
    # sp = State_Preparator(psi, qc)
    # time_evol = TimeEvolution(qc, psi, b, a)
    # meas = Measurement(qc)
    # qs = QuantumSimulator(qc, shots, psi)
    print("State preparation complete.")

    # Gaussian Potential Parameters
    V0 = 20.0
    sigma = 1.8
    x0 = len(psi) // 2
    # gp = GaussianPotential(V0, sigma, x0, psi, qc)
    # U = gp.create_potential_matrix()
    # Get the potential values directly from the GaussianPotential object
    x = np.arange(len(psi))
    dt=0.01
    # potential_values = [gp._potential(i) for i in x]
    # plotter = WavePlotter(num_sites=len(psi), potential_values=potential_values)

    plt.title("Gaussian Potential Well")
    plt.xlabel("Lattice Position")
    plt.ylabel("Potential U(x)")
    plt.grid(True)
    # print("\nPotential Matrix:\n", U)
    # print("Gaussian Potential instance created.")

    for t in range(time_steps):
        qc = QuantumCircuit(len(psi), len(psi))
        sp = State_Preparator(psi, qc)
        qc = sp.prepare()  
        meas = Measurement(qc)
        qs = QuantumSimulator(qc, shots, psi)
        # qc = gp.apply_potential(0.5*dt)                
        # qc = time_evol.apply_time_evolution(t)    # evolve t steps
        # qc = gp.apply_potential(0.5*dt)     # optional potential
        qc = meas.measure_all()             # add measurement gates

        counts = qs.run()                   # run simulator
        states_nonzero, amps_nonzero = qs.post_process(counts)

        print(f"\nAmps: {amps_nonzero}\n")
        print(f"\nProbs: {np.abs(amps_nonzero)**2}\n")

        # plotter.update(states_nonzero, amps_nonzero) 


main()