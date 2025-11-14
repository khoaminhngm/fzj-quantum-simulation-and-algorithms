import matplotlib.pyplot as plt
import numpy as np

class WavePlotter:
    def __init__(self, num_sites, potential_values=None):
        self.num_sites = num_sites
        self.potential = potential_values

        self.fig, self.ax = plt.subplots()
        self.x = np.arange(num_sites)
        
        # Bars for wave amplitude
        self.wave_bar = self.ax.bar(self.x, np.zeros(num_sites), color='blue', alpha=0.6)

        # Potential overlay (scaled to fit visually)
        if self.potential is not None:
            scaled = (self.potential - np.min(self.potential))
            scaled /= np.max(scaled)
            scaled *= 1.0   # scale relative to amplitudes
            self.ax.plot(self.x, scaled, color='red', linewidth=2, label="Potential")

        self.ax.set_xlabel("Position")
        self.ax.set_ylabel("Amplitude / Potential (scaled)")
        self.ax.set_title("Wave Function Evolution with Potential")
        self.ax.legend()
        plt.ion()
        plt.show()

    def update(self, states, amps):
        # Probability distribution (one value per lattice site)
        wave_prob = np.zeros(self.num_sites)

        for state, amp in zip(states, amps):
            if '1' not in state:
                continue
            pos = state.index('1')
            wave_prob[pos] = np.abs(amp)**2   # <-- probability

        # Update bar heights with probability values
        for bar, h in zip(self.wave_bar, wave_prob):
            bar.set_height(h)

        # Autoscale y-axis for visibility
        max_val = max(1e-10, np.max(wave_prob))
        self.ax.set_ylim(0, 1.2 * max_val)

        self.fig.canvas.draw()
        self.fig.canvas.flush_events()