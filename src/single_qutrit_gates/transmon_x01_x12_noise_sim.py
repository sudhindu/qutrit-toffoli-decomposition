"""
transmon_x01_x12_noise_sim.py

Noise simulation of the X01 and X12 single-qutrit gates on a driven
3-level transmon, quantifying how gate fidelity degrades under energy
relaxation and dephasing, via the Lindblad master equation in QuTiP.

Physical model
--------------
Unitary part: identical to the ideal (noiseless) simulation,
    H(t) = H0 + Omega(t) cos(wd t) * ( |0><1| + |1><0|
                                        + sqrt(2)(|1><2| + |2><1|) )
with the same Gaussian pulse calibrated to a pi-pulse area, for the
same representative transmon parameters (w01/2pi = 5.000 GHz,
alpha/2pi = 300 MHz).

Dissipative part: standard ladder-type relaxation,
    |1> -> |0>   at rate       Gamma1 = 1/T1
    |2> -> |1>   at rate     2*Gamma1
the factor-of-2 enhancement of the 2->1 decay rate is the *same*
sqrt(2)-dipole-matrix-element scaling (squared) already used for the
coherent drive coupling in Hd(t) -- i.e. it is the spontaneous-emission
counterpart of the same physical coupling, not an independent
assumption. Pure dephasing is included at rate Gamma_phi = Gamma1
(i.e. T_phi = T1), a standard default in the absence of an
independently measured T2, so that the single swept parameter T1
controls the entire noise budget.

Figure of merit
----------------
We use the "truth-table fidelity" F_TT: the population found in the
CORRECT target state, averaged over the 3 relevant computational
inputs, exactly the same metric already validated in the noiseless
simulation (transmon_x01_x12_sim.py) and the same style of metric used
by Galda et al. for a benchmarked CCNOT decomposition. As T1 ->
infinity, F_TT must converge exactly to the noiseless population-
transfer fidelities already reported (99.88% / 99.99% averaged over
the 3 inputs), which is used below as a built-in self-consistency
check on the code -- this was verified explicitly during development
and is why this metric is used here in preference to QuTiP's
higher-level Choi/process-fidelity routines, which were found to
disagree with direct state-population calculations for this specific
time-dependent, non-rotating-frame Hamiltonian.

Author: generated for Sudhindu Bikash Mandal, CQuERE.
Requires: numpy, qutip, matplotlib
"""

import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
import matplotlib as mpl

# ----------------------------------------------------------------------
# 1. Physical parameters -- identical to transmon_x01_x12_sim.py
# ----------------------------------------------------------------------
f01, alpha = 5.000, 0.300         # GHz
f12 = f01 - alpha

w01 = 2 * np.pi * f01
alpha_ang = 2 * np.pi * alpha
w12 = 2 * np.pi * f12

T_PULSE = 40.0            # ns
N_STEPS = 40000           # dt = 1 ps: resolves the 5 GHz carrier with ~200 pts/cycle
tlist = np.linspace(0.0, T_PULSE, N_STEPS)

g0, g1, g2 = qt.basis(3, 0), qt.basis(3, 1), qt.basis(3, 2)
proj1, proj2 = g1 * g1.dag(), g2 * g2.dag()
b01 = g0 * g1.dag() + g1 * g0.dag()
b12 = g1 * g2.dag() + g2 * g1.dag()

H0 = w01 * proj1 + (2 * w01 - alpha_ang) * proj2
Hd_op = b01 + np.sqrt(2) * b12

def gaussian_shape(t, t0, sigma):
    return np.exp(-0.5 * ((t - t0) / sigma) ** 2)

def build_envelope(area_target, coupling_factor=1.0):
    t0, sigma = T_PULSE / 2, T_PULSE / 6.0
    shape = gaussian_shape(tlist, t0, sigma)
    raw_area = np.trapezoid(shape, tlist)
    return (area_target / (coupling_factor * raw_area)) * shape

Omega_X01 = build_envelope(np.pi, 1.0)
Omega_X12 = build_envelope(np.pi, np.sqrt(2))

def make_H(Omega_array, wd):
    return qt.QobjEvo([H0, [Hd_op, Omega_array * np.cos(wd * tlist)]], tlist=tlist)

H_X01 = make_H(Omega_X01, w01)
H_X12 = make_H(Omega_X12, w12)

# Integrator options: a coarse 2-point output request ([0, T_PULSE]) still
# requires an explicit step-size cap so the adaptive solver resolves the
# ~5 GHz carrier oscillation within a single internal integration call.
SOLVER_OPTS = {"nsteps": 200000, "max_step": 0.01}

# ----------------------------------------------------------------------
# 2. Collapse operators for a given T1 (in ns)
# ----------------------------------------------------------------------
def collapse_ops(T1_ns):
    Gamma1 = 1.0 / T1_ns
    c_relax_10 = np.sqrt(Gamma1) * (g0 * g1.dag())        # |1> -> |0>
    c_relax_21 = np.sqrt(2 * Gamma1) * (g1 * g2.dag())    # |2> -> |1>, enhanced x2
    dephasing_op = qt.Qobj(np.diag([0.0, 1.0, 2.0]))
    c_deph = np.sqrt(2 * Gamma1) * dephasing_op           # T_phi = T1
    return [c_relax_10, c_relax_21, c_deph]

# ----------------------------------------------------------------------
# 3. Truth-table fidelity: population in the correct target state,
#    averaged over the 3 relevant computational-basis inputs.
# ----------------------------------------------------------------------
basis_kets = [g0, g1, g2]
targets_X01 = [g1, g0, g2]   # X01: |0>->|1>, |1>->|0>, |2>->|2>
targets_X12 = [g0, g2, g1]   # X12: |0>->|0>, |1>->|2>, |2>->|1>

def truth_table_fidelity(H, targets, c_ops):
    fids = []
    for psi0, target in zip(basis_kets, targets):
        rho0 = psi0 * psi0.dag()
        res = qt.mesolve(H, rho0, [0, T_PULSE], c_ops=c_ops, options=SOLVER_OPTS)
        rho_out = res.states[-1]
        fids.append(qt.expect(target * target.dag(), rho_out))
    return float(np.real(np.mean(fids)))

# ----------------------------------------------------------------------
# 4. Self-consistency check: the T1 -> infinity (zero-noise) limit must
#    reproduce the noiseless population-transfer fidelities already
#    reported in transmon_x01_x12_sim.py (99.88% / 99.99% averaged).
# ----------------------------------------------------------------------
f_check_X01 = truth_table_fidelity(H_X01, targets_X01, [])
f_check_X12 = truth_table_fidelity(H_X12, targets_X12, [])
print(f"Zero-noise check: F_TT(X01) = {f_check_X01:.6f}  (expect ~0.998817)")
print(f"Zero-noise check: F_TT(X12) = {f_check_X12:.6f}  (expect ~0.999933)")

# ----------------------------------------------------------------------
# 5. Sweep T1 and compute the truth-table fidelity for both gates
# ----------------------------------------------------------------------
T1_us_list = np.logspace(0, 3, 16)   # 1 us to 1000 us
T1_ns_list = T1_us_list * 1e3

fid_X01, fid_X12 = [], []
for T1_ns in T1_ns_list:
    c_ops = collapse_ops(T1_ns)
    f01_ = truth_table_fidelity(H_X01, targets_X01, c_ops)
    f12_ = truth_table_fidelity(H_X12, targets_X12, c_ops)
    fid_X01.append(f01_)
    fid_X12.append(f12_)
    print(f"T1 = {T1_ns/1e3:8.3f} us   F_TT(X01) = {f01_:.6f}   F_TT(X12) = {f12_:.6f}")

fid_X01 = np.array(fid_X01)
fid_X12 = np.array(fid_X12)
infid_X01 = 1 - fid_X01
infid_X12 = 1 - fid_X12

print(f"\nSanity check at largest T1 = {T1_us_list[-1]:.0f} us: "
      f"F_TT(X01) = {fid_X01[-1]:.6f} (noiseless floor {f_check_X01:.6f}), "
      f"F_TT(X12) = {fid_X12[-1]:.6f} (noiseless floor {f_check_X12:.6f}).")

# ----------------------------------------------------------------------
# 6. Plot
# ----------------------------------------------------------------------
mpl.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "axes.linewidth": 1.0,
})

NAVY = "#0B1F3A"
RED = "#C0442C"

fig, ax = plt.subplots(figsize=(6.0, 4.6))

ax.plot(T1_us_list, infid_X01, "o-", color=NAVY, lw=1.8, ms=6,
        label=r"$X_{01}$ gate")
ax.plot(T1_us_list, infid_X12, "s-", color=RED, lw=1.8, ms=6,
        label=r"$X_{12}$ gate")

ax.axhline(1 - f_check_X01, color=NAVY, lw=1.0, ls=":", alpha=0.7,
           label=r"noiseless floor ($X_{01}$)")
ax.axhline(1 - f_check_X12, color=RED, lw=1.0, ls=":", alpha=0.7,
           label=r"noiseless floor ($X_{12}$)")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$T_1$ ($\mu$s)")
ax.set_ylabel(r"Gate infidelity, $1 - F_{TT}$")
ax.set_title(r"Noise-induced infidelity of $X_{01}$ and $X_{12}$ vs. $T_1$"
             "\n" r"(pure dephasing fixed at $T_\phi = T_1$)", fontsize=12)
ax.legend(loc="lower left", fontsize=9.5, frameon=True, ncol=1)
ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.6)

fig.tight_layout()
fig.savefig("X01_X12_noise_fidelity.png", dpi=220)
print("\nX01_X12_noise_fidelity.png")

# ----------------------------------------------------------------------
# 7. Save numeric data for independent verification / re-plotting
# ----------------------------------------------------------------------
np.savetxt("noise_fidelity_data.csv",
           np.column_stack([T1_us_list, fid_X01, fid_X12, infid_X01, infid_X12]),
           header="T1_us, F_TT_X01, F_TT_X12, infid_X01, infid_X12", delimiter=",")
print("noise_fidelity_data.csv")
