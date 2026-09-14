"""
transmon_x01_x12_sim.py

Full lab-frame (no rotating-wave pre-approximation) simulation of a
driven 3-level transmon, used to numerically validate the physical
realization of the X01 and X12 single-qutrit gates described in the
paper:

    H0 = w01 |1><1| + (2*w01 - alpha) |2><2|

    Hd(t) = Omega(t) * cos(wd*t) * ( |0><1| + |1><0|
                                      + sqrt(2) * (|1><2| + |2><1|) )

    H(t)  = H0 + Hd(t)

For the X01 gate, wd = w01 and the pulse area satisfies
    integral( Omega(t) dt ) = pi
For the X12 gate, wd = w12 = w01 - alpha and
    integral( sqrt(2) * Omega(t) dt ) = pi

The script does NOT invoke the rotating-wave approximation anywhere:
the Schrodinger equation is integrated in the lab frame, with the
carrier cos(wd t) included explicitly. This means any leakage caused
by the counter-rotating / off-resonant terms that the RWA discards
analytically is captured numerically. The purpose of the simulation
is precisely to confirm that this leakage is negligible for realistic
parameters (alpha >> Omega_max), i.e. to validate --not assume-- the
analytic isolation proof given in the paper.

Author: generated for Sudhindu Bikash Mandal, CQuERE.
Requires: numpy, qutip, matplotlib
"""

import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
import matplotlib as mpl

# ----------------------------------------------------------------------
# 1. Physical parameters (transmon-like, in linear frequency units GHz;
#    time is in ns, so that frequency [GHz] * time [ns] is dimensionless
#    once multiplied by 2*pi to convert to angular frequency).
# ----------------------------------------------------------------------
f01    = 5.000       # 0<->1 transition frequency,   GHz
alpha  = 0.300        # anharmonicity alpha = w01 - w12,  GHz  (>0 by convention used in the paper)
f12    = f01 - alpha  # 1<->2 transition frequency,   GHz

w01 = 2 * np.pi * f01      # angular frequency, rad/ns
w12 = 2 * np.pi * f12      # angular frequency, rad/ns
alpha_ang = 2 * np.pi * alpha

T_PULSE = 40.0    # total pulse duration, ns
N_STEPS = 40000   # time samples (dt = 1 ps; carrier period at 5 GHz is 200 ps, so this
                  # resolves the fastest oscillation with ~200 samples/cycle)

tlist = np.linspace(0.0, T_PULSE, N_STEPS)
dt = tlist[1] - tlist[0]

# ----------------------------------------------------------------------
# 2. Operators (3-level Hilbert space: |0>, |1>, |2>)
# ----------------------------------------------------------------------
g0, g1, g2 = qt.basis(3, 0), qt.basis(3, 1), qt.basis(3, 2)

proj0 = g0 * g0.dag()
proj1 = g1 * g1.dag()
proj2 = g2 * g2.dag()

b01 = g0 * g1.dag() + g1 * g0.dag()   # |0><1| + |1><0|
b12 = g1 * g2.dag() + g2 * g1.dag()   # |1><2| + |2><1|

H0 = w01 * proj1 + (2 * w01 - alpha_ang) * proj2
Hd_op = b01 + np.sqrt(2) * b12         # the operator multiplying Omega(t)*cos(wd t)

# ----------------------------------------------------------------------
# 3. Pulse envelope: truncated Gaussian, a standard experimentally
#    realistic shape (smooth turn-on/off avoids spectral leakage).
#    The amplitude is solved numerically so that the *exact* discretized
#    integral gives the requested pulse area (pi, or pi/sqrt(2)).
# ----------------------------------------------------------------------
def gaussian_shape(t, t0, sigma):
    """Unit-height Gaussian envelope centered at t0."""
    return np.exp(-0.5 * ((t - t0) / sigma) ** 2)

def build_envelope(area_target, coupling_factor=1.0):
    """
    Returns Omega(t) sampled on `tlist` such that
        coupling_factor * integral(Omega(t) dt) = area_target
    using a truncated Gaussian shape centered in the pulse window.
    """
    t0 = T_PULSE / 2
    sigma = T_PULSE / 6.0   # +/-3 sigma fits inside [0, T_PULSE]
    shape = gaussian_shape(tlist, t0, sigma)
    raw_area = np.trapezoid(shape, tlist)
    amplitude = area_target / (coupling_factor * raw_area)
    return amplitude * shape

# X01 gate: drive at w01, area(Omega) = pi
Omega_X01 = build_envelope(area_target=np.pi, coupling_factor=1.0)

# X12 gate: drive at w12, area(sqrt(2)*Omega) = pi  =>  area(Omega) = pi/sqrt(2)
Omega_X12 = build_envelope(area_target=np.pi, coupling_factor=np.sqrt(2))

# ----------------------------------------------------------------------
# 4. Time-dependent Hamiltonian builder (array-based QobjEvo -> fast)
# ----------------------------------------------------------------------
def make_hamiltonian(Omega_array, wd):
    drive_coeff = Omega_array * np.cos(wd * tlist)
    H = [H0, [Hd_op, drive_coeff]]
    return qt.QobjEvo(H, tlist=tlist)

# ----------------------------------------------------------------------
# 5. Simulate the gate starting from each computational basis state and
#    record the full population dynamics P0(t), P1(t), P2(t).
# ----------------------------------------------------------------------
def simulate_gate(Omega_array, wd):
    H = make_hamiltonian(Omega_array, wd)
    results = {}
    for label, psi0 in zip(["0", "1", "2"], [g0, g1, g2]):
        res = qt.sesolve(H, psi0, tlist, e_ops=[proj0, proj1, proj2])
        results[label] = {
            "P0": res.expect[0],
            "P1": res.expect[1],
            "P2": res.expect[2],
        }
    return results

# ----------------------------------------------------------------------
# 6. Build the final-time population-transfer ("truth") table:
#    T[i][j] = P(measure j | start in i), at t = T_PULSE
# ----------------------------------------------------------------------
def transfer_table(results):
    table = np.zeros((3, 3))
    for i, label in enumerate(["0", "1", "2"]):
        table[i, 0] = results[label]["P0"][-1]
        table[i, 1] = results[label]["P1"][-1]
        table[i, 2] = results[label]["P2"][-1]
    return table

print("Running X01 gate simulation (drive resonant with |0><->|1>) ...")
res_X01 = simulate_gate(Omega_X01, w01)
table_X01 = transfer_table(res_X01)

print("Running X12 gate simulation (drive resonant with |1><->|2>) ...")
res_X12 = simulate_gate(Omega_X12, w12)
table_X12 = transfer_table(res_X12)

np.set_printoptions(precision=4, suppress=True)
print("\nX01 population-transfer table  (rows = input, cols = output |0>,|1>,|2>):")
print(table_X01)
print("\nX12 population-transfer table  (rows = input, cols = output |0>,|1>,|2>):")
print(table_X12)

# ----------------------------------------------------------------------
# 7. Plotting
# ----------------------------------------------------------------------
mpl.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.linewidth": 0.9,
})

NAVY = "#0B1F3A"
BLUE = "#3B82C4"
RED = "#C0442C"
GREEN = "#3E8E5A"
GRAY = "#5A6B7D"

def plot_gate(Omega_array, wd, results, table, gate_name, spectator_label, fname):
    fig, axes = plt.subplots(4, 1, figsize=(6.4, 9.2), sharex=True,
                              gridspec_kw={"height_ratios": [1, 1.3, 1.3, 1.3]})

    # --- Panel A: pulse envelope ---
    ax = axes[0]
    ax.plot(tlist, Omega_array / (2 * np.pi) * 1000, color=NAVY, lw=1.8)
    ax.set_ylabel(r"$\Omega(t)/2\pi$ (MHz)")
    ax.set_title(f"{gate_name} gate drive envelope  "
                 f"($\\omega_d/2\\pi$ = {wd/(2*np.pi):.3f} GHz)", fontsize=12)
    ax.axhline(0, color="black", lw=0.6)

    labels = ["0", "1", "2"]
    colors = {"P0": BLUE, "P1": RED, "P2": GREEN}
    for k, init_label in enumerate(labels):
        ax2 = axes[k + 1]
        r = results[init_label]
        ax2.plot(tlist, r["P0"], color=colors["P0"], lw=1.6, label=r"$P_{|0\rangle}(t)$")
        ax2.plot(tlist, r["P1"], color=colors["P1"], lw=1.6, label=r"$P_{|1\rangle}(t)$")
        ax2.plot(tlist, r["P2"], color=colors["P2"], lw=1.6, label=r"$P_{|2\rangle}(t)$")
        ax2.set_ylim(-0.05, 1.05)
        ax2.set_ylabel("Population")
        tag = " (spectator)" if init_label == spectator_label else ""
        ax2.text(0.02, 0.92, f"Input: $|{init_label}\\rangle$" + tag,
                  transform=ax2.transAxes, fontsize=10.5, va="top",
                  color=NAVY, fontweight="bold")
        final_str = (f"Final: $P_0$={r['P0'][-1]:.3f}, "
                      f"$P_1$={r['P1'][-1]:.3f}, $P_2$={r['P2'][-1]:.3f}")
        ax2.text(0.02, 0.08, final_str, transform=ax2.transAxes,
                  fontsize=9.5, va="bottom", color=GRAY)
        if k == 0:
            ax2.legend(loc="upper right", fontsize=9, ncol=3, frameon=False)

    axes[-1].set_xlabel("Time (ns)")
    fig.tight_layout()
    fig.savefig(fname, dpi=220)
    plt.close(fig)
    print(f"Saved {fname}")

plot_gate(Omega_X01, w01, res_X01, table_X01,
          gate_name=r"$X_{01}$", spectator_label="2",
          fname="X01_gate_simulation.png")

plot_gate(Omega_X12, w12, res_X12, table_X12,
          gate_name=r"$X_{12}$", spectator_label="0",
          fname="X12_gate_simulation.png")

# ----------------------------------------------------------------------
# 8. Save numeric tables to a text file for direct inclusion in the
#    paper's Supplemental Material.
# ----------------------------------------------------------------------
with open("transfer_tables.txt", "w") as f:
    f.write("X01 gate population-transfer table (rows=input, cols=output |0>,|1>,|2>)\n")
    f.write(np.array2string(table_X01, precision=5) + "\n\n")
    f.write("X12 gate population-transfer table (rows=input, cols=output |0>,|1>,|2>)\n")
    f.write(np.array2string(table_X12, precision=5) + "\n")

print("\nDone. Tables saved to transfer_tables.txt")
