# Deutsch's Algorithm

Determining if \( f(0) = f(1) \) with one run on two qubits.

**Start with:**
$$
|\psi_0\rangle = |0\rangle|1\rangle
$$

**Apply Hadamards:**
$$
|\psi_1\rangle = (H\otimes H)|01\rangle 
= \tfrac{1}{2}(|00\rangle - |01\rangle + |10\rangle - |11\rangle)
$$

**Oracle action:**
$$
U_f|x,y\rangle = |x, y\oplus f(x)\rangle
$$
$$
U_f|\psi_1\rangle = \tfrac{1}{2}\big((-1)^{f(0)}|0\rangle + (-1)^{f(1)}|1\rangle\big)(|0\rangle - |1\rangle)
$$

**First-qubit state:**
$$
|\psi_x\rangle = \tfrac{1}{\sqrt{2}}\big((-1)^{f(0)}|0\rangle + (-1)^{f(1)}|1\rangle\big)
$$

**Apply Hadamard:**
$$
H|\psi_x\rangle = 
\tfrac{1}{2}\big([(-1)^{f(0)}+(-1)^{f(1)}]|0\rangle + [(-1)^{f(0)}-(-1)^{f(1)}]|1\rangle\big)
$$

**Measure first qubit:**
$$
\boxed{
\begin{cases}
|0\rangle & f(0)=f(1)\ \text{(constant)}\\[4pt]
|1\rangle & f(0)\neq f(1)\ \text{(balanced)}
\end{cases}
}
$$
