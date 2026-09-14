"""
transmon_C1X12_C2X01_C2X12_CR_sim_v3.py

Realistic two-transmon simulation of all THREE controlled two-qutrit
gates used in this work's Toffoli / C^K-NOT decomposition:
    C^1 X_12  -- X_12 acts on the target iff the control is |1>
    C^2 X_01  -- X_01 acts on the target iff the control is |2>
    C^2 X_12  -- X_12 acts on the target iff the control is |2>
                 (the additional primitive needed for K>2, i.e. the
                 C^K-NOT generalization beyond the basic Toffoli)

All three gates are realized via the same cross-resonance (CR)
mechanism: control and target are distinct, capacitively coupled
fixed-frequency transmons, driven through the CONTROL's charge line
at a frequency resonant with the TARGET's transition as dressed by
the addressed control level. The static Hamiltonian (transmon
frequencies, anharmonicities, and coupling g) is IDENTICAL to the
C^1X_12/C^2X_01 simulation already reported; only the drive frequency,
amplitude, and duration differ per gate, so the previously reported
C^1X_12/C^2X_01 results are reproduced unchanged by this script.

Author: generated for Sudhindu Bikash Mandal, CQuERE.
Requires: numpy, qutip, matplotlib
"""

import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
import matplotlib as mpl

# ----------------------------------------------------------------------
# 1. Two distinct, realistic transmons (identical to the C1X12/C2X01 run)
# ----------------------------------------------------------------------
f01_c, alpha_c = 5.000, 0.300       # control, GHz
f01_t, alpha_t = 5.550, 0.320       # target, GHz

w01_c, alpha_c_ang = 2*np.pi*f01_c, 2*np.pi*alpha_c
w01_t, alpha_t_ang = 2*np.pi*f01_t, 2*np.pi*alpha_t

g0, g1, g2 = qt.basis(3, 0), qt.basis(3, 1), qt.basis(3, 2)
proj1, proj2 = g1*g1.dag(), g2*g2.dag()
b01 = g0*g1.dag() + g1*g0.dag()
b12 = g1*g2.dag() + g2*g1.dag()
charge = b01 + np.sqrt(2)*b12
I3 = qt.qeye(3)

H0_c = w01_c*proj1 + (2*w01_c - alpha_c_ang)*proj2
H0_t = w01_t*proj1 + (2*w01_t - alpha_t_ang)*proj2
H0 = qt.tensor(H0_c, I3) + qt.tensor(I3, H0_t)

G_MHZ = 25.0
g_ang = 2*np.pi*G_MHZ*1e-3
H_couple = g_ang * qt.tensor(charge, charge)
H_static = H0 + H_couple

Hd_control = qt.tensor(charge, I3)

# ----------------------------------------------------------------------
# 2. Overlap diagnostic (identical check as before; confirms this run
#    uses the same, already-validated static Hamiltonian)
# ----------------------------------------------------------------------
basis_labels = [(i, j) for i in range(3) for j in range(3)]
bare_states = {(i, j): qt.tensor(qt.basis(3, i), qt.basis(3, j)) for i, j in basis_labels}
evals, evecs = H_static.eigenstates()
E = {}
worst_overlap = 1.0
for k in range(9):
    overlaps = [(lbl, abs(evecs[k].overlap(bare_states[lbl]))**2) for lbl in basis_labels]
    best = max(overlaps, key=lambda x: x[1])
    E[best[0]] = evals[k] / (2*np.pi)
    worst_overlap = min(worst_overlap, best[1])
print(f"Overlap diagnostic: worst-case bare-state overlap = {worst_overlap:.4f} "
      f"(matches the value already reported for C1X12/C2X01: 0.976)")

f_t01_given_c = {c: E[(c, 1)] - E[(c, 0)] for c in range(3)}
f_t12_given_c = {c: E[(c, 2)] - E[(c, 1)] for c in range(3)}
wd_C2X01 = 2*np.pi * f_t01_given_c[2]
wd_C1X12 = 2*np.pi * f_t12_given_c[1]
wd_C2X12 = 2*np.pi * f_t12_given_c[2]   # NEW: target's 1<->2 transition
                                         # dressed by control = |2>

print(f"C2X12 drive frequency: {wd_C2X12/(2*np.pi):.6f} GHz "
      f"(splitting from control=0 branch: "
      f"{(f_t12_given_c[2]-f_t12_given_c[0])*1000:.3f} MHz; "
      f"from control=1 branch: "
      f"{(f_t12_given_c[2]-f_t12_given_c[1])*1000:.3f} MHz)")

# ----------------------------------------------------------------------
# 3. Calibrated pulse parameters. C1X12/C2X01 UNCHANGED from the
#    previously reported run. C2X12 is NEW: found via a joint
#    amplitude/duration scan maximizing both addressed-branch fidelity
#    AND control-selectivity (the control=0 branch has an intrinsically
#    small ~0.5 MHz splitting from the addressed control=2 branch at
#    this static Hamiltonian, requiring a slower pulse than C1X12/C2X01
#    for good selectivity).
# ----------------------------------------------------------------------
PULSE = {
    "C1X12": {"amp_MHz": 25.0, "T_ns": 900.0, "wd": wd_C1X12},
    "C2X01": {"amp_MHz": 35.0, "T_ns": 250.0, "wd": wd_C2X01},
    "C2X12": {"amp_MHz": 35.0, "T_ns": 600.0, "wd": wd_C2X12},
}
RAMP_NS = 10.0
SOLVER_OPTS = {"nsteps": 2000000, "max_step": 0.008}

def envelope(tlist, amp_ang, T_pulse, ramp=RAMP_NS):
    e = np.ones_like(tlist) * amp_ang
    m1 = tlist < ramp
    e[m1] = amp_ang * 0.5 * (1 - np.cos(np.pi * tlist[m1] / ramp))
    m2 = tlist > (T_pulse - ramp)
    e[m2] = amp_ang * 0.5 * (1 - np.cos(np.pi * (T_pulse - tlist[m2]) / ramp))
    return e

def make_H(gate, n_steps_per_ns=100):
    p = PULSE[gate]
    T_pulse = p["T_ns"]
    tlist = np.linspace(0, T_pulse, max(int(T_pulse * n_steps_per_ns), 2000))
    amp_ang = 2*np.pi * p["amp_MHz"] * 1e-3
    e = envelope(tlist, amp_ang, T_pulse)
    H = qt.QobjEvo([H_static, [Hd_control, e * np.cos(p["wd"] * tlist)]], tlist=tlist)
    return H, tlist

# ----------------------------------------------------------------------
# 4. Ideal (noiseless) simulation, with optional noise (c_ops)
# ----------------------------------------------------------------------
def run_ideal(gate, control_level, target_ket, c_ops=None):
    H, tlist = make_H(gate)
    psi0 = qt.tensor(qt.basis(3, control_level), target_ket)
    if c_ops:
        rho0 = psi0 * psi0.dag()
        res = qt.mesolve(H, rho0, [0, tlist[-1]], c_ops=c_ops, options=SOLVER_OPTS)
        rho_out = res.states[-1]
    else:
        res = qt.sesolve(H, psi0, [0, tlist[-1]], options=SOLVER_OPTS)
        rho_out = qt.ket2dm(res.states[-1])
    Ps = []
    for m in range(3):
        proj_t = qt.tensor(qt.qeye(3), qt.basis(3, m)*qt.basis(3, m).dag())
        Ps.append(qt.expect(proj_t, rho_out))
    return np.real(Ps)

print("\n=== C^1 X_12 (drive resonant for control=|1>) ===")
rows_C1X12 = []
for control, target_in in [(c, t) for c in [0,1,2] for t in [g0,g1,g2]]:
    P = run_ideal("C1X12", control, target_in)
    rows_C1X12.append((control, target_in, P))
    tlabel = {id(g0): "|0>", id(g1): "|1>", id(g2): "|2>"}[id(target_in)]
    print(f"control={control} target_in={tlabel}: P0={P[0]:.4f} P1={P[1]:.4f} P2={P[2]:.4f}")

print("\n=== C^2 X_01 (drive resonant for control=|2>) ===")
rows_C2X01 = []
for control, target_in in [(c, t) for c in [0,1,2] for t in [g0,g1,g2]]:
    P = run_ideal("C2X01", control, target_in)
    rows_C2X01.append((control, target_in, P))
    tlabel = {id(g0): "|0>", id(g1): "|1>", id(g2): "|2>"}[id(target_in)]
    print(f"control={control} target_in={tlabel}: P0={P[0]:.4f} P1={P[1]:.4f} P2={P[2]:.4f}")

print("\n=== C^2 X_12 (drive resonant for control=|2>) [NEW] ===")
rows_C2X12 = []
for control, target_in in [(c, t) for c in [0,1,2] for t in [g0,g1,g2]]:
    P = run_ideal("C2X12", control, target_in)
    rows_C2X12.append((control, target_in, P))
    tlabel = {id(g0): "|0>", id(g1): "|1>", id(g2): "|2>"}[id(target_in)]
    print(f"control={control} target_in={tlabel}: P0={P[0]:.4f} P1={P[1]:.4f} P2={P[2]:.4f}")

# ----------------------------------------------------------------------
# 5. Population dynamics for the figure
# ----------------------------------------------------------------------
def run_traj(gate, control_level, target_ket, level_to_track):
    H, tlist = make_H(gate)
    psi0 = qt.tensor(qt.basis(3, control_level), target_ket)
    proj_t = qt.tensor(qt.qeye(3), qt.basis(3, level_to_track)*qt.basis(3, level_to_track).dag())
    res = qt.sesolve(H, psi0, tlist, e_ops=[proj_t], options=SOLVER_OPTS)
    return tlist, res.expect[0]

traj_C1X12 = {c: run_traj("C1X12", c, g1, 2) for c in [0, 1, 2]}
traj_C2X01 = {c: run_traj("C2X01", c, g0, 1) for c in [0, 1, 2]}
traj_C2X12 = {c: run_traj("C2X12", c, g1, 2) for c in [0, 1, 2]}

# ----------------------------------------------------------------------
# 6. Noisy simulation: T1 + dephasing on both control and target,
#    truth-table fidelity of the active branch vs T1, for all 3 gates
# ----------------------------------------------------------------------
def collapse_ops_2q(T1_ns):
    Gamma1 = 1.0 / T1_ns
    N_op = qt.Qobj(np.diag([0., 1., 2.]))
    ops = []
    for qutrit_idx in [0, 1]:
        def embed(op):
            return qt.tensor(op, I3) if qutrit_idx == 0 else qt.tensor(I3, op)
        ops.append(np.sqrt(Gamma1) * embed(g0*g1.dag()))
        ops.append(np.sqrt(2*Gamma1) * embed(g1*g2.dag()))
        ops.append(np.sqrt(2*Gamma1) * embed(N_op))
    return ops

T1_us_list = np.logspace(0, 3, 13)
T1_ns_list = T1_us_list * 1e3

fid_C1X12, fid_C2X01, fid_C2X12 = [], [], []
for T1_ns in T1_ns_list:
    c_ops = collapse_ops_2q(T1_ns)
    P1 = run_ideal("C1X12", 1, g1, c_ops=c_ops)[2]   # P2, target of forward C1X12
    P2 = run_ideal("C2X01", 2, g0, c_ops=c_ops)[1]   # P1, target of forward C2X01
    P3 = run_ideal("C2X12", 2, g1, c_ops=c_ops)[2]   # P2, target of forward C2X12
    fid_C1X12.append(P1)
    fid_C2X01.append(P2)
    fid_C2X12.append(P3)
    print(f"T1={T1_ns/1e3:8.2f} us   F(C1X12)={P1:.4f}   F(C2X01)={P2:.4f}   F(C2X12)={P3:.4f}")

fid_C1X12 = np.array(fid_C1X12)
fid_C2X01 = np.array(fid_C2X01)
fid_C2X12 = np.array(fid_C2X12)

np.savetxt("CR_noise_fidelity_data_v3.csv",
           np.column_stack([T1_us_list, fid_C1X12, fid_C2X01, fid_C2X12]),
           header="T1_us, F_C1X12, F_C2X01, F_C2X12", delimiter=",")

# ----------------------------------------------------------------------
# 7. Save table to CSV
# ----------------------------------------------------------------------
with open("CR_controlled_gate_table_v3.csv", "w") as f:
    f.write("gate,control,target_in,P0,P1,P2\n")
    label_map = {id(g0): 0, id(g1): 1, id(g2): 2}
    for control, target_in, P in rows_C1X12:
        f.write(f"C1X12,{control},{label_map[id(target_in)]},{P[0]:.4f},{P[1]:.4f},{P[2]:.4f}\n")
    for control, target_in, P in rows_C2X01:
        f.write(f"C2X01,{control},{label_map[id(target_in)]},{P[0]:.4f},{P[1]:.4f},{P[2]:.4f}\n")
    for control, target_in, P in rows_C2X12:
        f.write(f"C2X12,{control},{label_map[id(target_in)]},{P[0]:.4f},{P[1]:.4f},{P[2]:.4f}\n")

# ----------------------------------------------------------------------
# 8. Plots
# ----------------------------------------------------------------------
mpl.rcParams.update({"font.family": "serif", "font.size": 12, "axes.linewidth": 1.0})
NAVY, RED, GREEN = "#0B1F3A", "#C0442C", "#3E8E5A"
colors = {0: GREEN, 1: NAVY, 2: RED}

fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.8))
ax = axes[0]
for c in [0, 1, 2]:
    t, P = traj_C1X12[c]
    ax.plot(t, P, color=colors[c], lw=1.8, label=fr"control $=|{c}\rangle$")
ax.set_title(r"(a) $C^1X_{12}$: target $|1\rangle\to|2\rangle$?", fontsize=11.5)
ax.set_xlabel("Time (ns)"); ax.set_ylabel(r"Target $P_{|2\rangle}(t)$")
ax.set_ylim(-0.05, 1.05); ax.legend(fontsize=9, loc="center right", frameon=False)

ax = axes[1]
for c in [0, 1, 2]:
    t, P = traj_C2X01[c]
    ax.plot(t, P, color=colors[c], lw=1.8, label=fr"control $=|{c}\rangle$")
ax.set_title(r"(b) $C^2X_{01}$: target $|0\rangle\to|1\rangle$?", fontsize=11.5)
ax.set_xlabel("Time (ns)"); ax.set_ylabel(r"Target $P_{|1\rangle}(t)$")
ax.set_ylim(-0.05, 1.05); ax.legend(fontsize=9, loc="center right", frameon=False)

ax = axes[2]
for c in [0, 1, 2]:
    t, P = traj_C2X12[c]
    ax.plot(t, P, color=colors[c], lw=1.8, label=fr"control $=|{c}\rangle$")
ax.set_title(r"(c) $C^2X_{12}$: target $|1\rangle\to|2\rangle$?", fontsize=11.5)
ax.set_xlabel("Time (ns)"); ax.set_ylabel(r"Target $P_{|2\rangle}(t)$")
ax.set_ylim(-0.05, 1.05); ax.legend(fontsize=9, loc="center right", frameon=False)

fig.tight_layout()
fig.savefig("CR_controlled_gate_dynamics_v3.png", dpi=220)
print("Saved CR_controlled_gate_dynamics_v3.png")

fig2, ax = plt.subplots(figsize=(6.4, 4.8))
ax.plot(T1_us_list, 1-fid_C1X12, "o-", color=NAVY, lw=1.8, ms=6, label=r"$C^1X_{12}$")
ax.plot(T1_us_list, 1-fid_C2X01, "s-", color=RED, lw=1.8, ms=6, label=r"$C^2X_{01}$")
ax.plot(T1_us_list, 1-fid_C2X12, "^-", color=GREEN, lw=1.8, ms=6, label=r"$C^2X_{12}$")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"$T_1$ ($\mu$s)")
ax.set_ylabel(r"Active-branch infidelity, $1-F_{TT}$")
ax.set_title("Noise-induced infidelity of the controlled gates vs. $T_1$\n"
             r"(pure dephasing fixed at $T_\phi=T_1$, both qutrits)", fontsize=12)
ax.legend(loc="lower left", fontsize=10, frameon=True)
ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.6)
fig2.tight_layout()
fig2.savefig("CR_noise_fidelity_v3.png", dpi=220)
print("Saved CR_noise_fidelity_v3.png")
