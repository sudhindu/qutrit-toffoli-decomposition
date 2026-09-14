"""
toffoli_noise_comparison_v2.py

Projected Toffoli-gate fidelity vs. T1, updated to:
  (1) use the AUDITED physical pulse count for Galda et al. (12, not
      the entangling-pulse-only count of 4 originally reported in
      their Table III row), derived by explicitly decomposing their
      C^1X_12 circuit block (X_- . CNOT . X_+ = X_12, verified
      algebraically) into its constituent physical pulses:
          X_-  = 2 single-qutrit pulses
          CNOT = 1 two-qutrit pulse
          X_+  = 2 single-qutrit pulses
      -> 5 pulses per C^1X_12-equivalent block, x2 (fold + unfold)
         = 10, plus their middle C^2NOT (2 pulses, per their own
         stated "two consecutive CNOT gates" realization) = 12 total.
  (2) recompute the common per-pulse duration using this corrected
      denominator: t_pulse = 1.593 us / 12 = 132.75 ns (previously
      1.593 us / 4 = 398.25 ns, based on the uncorrected count).
  (3) add a new "hybrid" construction: our own native, self-inverse
      X_01/X_12 gates used to bracket Galda's existing CX primitive
      (in place of their composite, cyclic X_-/X_+), while leaving
      their already-efficient 2-pulse C^2NOT trick untouched. This
      isolates the pulse-count benefit attributable specifically to
      subspace-selective vs. cyclic single-qutrit gates, independent
      of adopting the fully-native CR two-qutrit mechanism.
        2 x (X_01/X_12 + CX + X_12/X_01) = 2 x 3 = 6 pulses
        + Galda's own middle C^2NOT (2 pulses, unchanged)
        = 8 total.

Requires: numpy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# ----------------------------------------------------------------------
# 1. Physical pulse counts (audited where noted)
# ----------------------------------------------------------------------
constructions = {
    "This work":              {"N_pulses": 3,  "color": "#0B1F3A", "marker": "o"},
    "This work + Galda CX":   {"N_pulses": 8,  "color": "#7A4FA3", "marker": "P"},
    "Gokhale et al.":         {"N_pulses": 6,  "color": "#3B82C4", "marker": "^"},
    "Galda et al. ": {"N_pulses": 12, "color": "#C0442C", "marker": "s"},
    "Nikolaeva et al.":       {"N_pulses": 4,  "color": "#3E8E5A", "marker": "D"},
}

T_PULSE_NS = 1593.0 / 12   # 132.75 ns; corrected common per-pulse duration

# ----------------------------------------------------------------------
# 2. Fidelity model (unchanged formula, corrected inputs)
# ----------------------------------------------------------------------
def fidelity(T1_ns, N_pulses, t_pulse_ns=T_PULSE_NS):
    per_pulse_infidelity = (2.0/3.0) * (t_pulse_ns / T1_ns)
    return np.clip(1.0 - per_pulse_infidelity, 0.0, 1.0) ** N_pulses

T1_us_list = np.logspace(0, 3, 200)
T1_ns_list = T1_us_list * 1e3

# ----------------------------------------------------------------------
# 3. Plot
# ----------------------------------------------------------------------
mpl.rcParams.update({"font.family": "serif", "font.size": 12, "axes.linewidth": 1.0})

fig, ax = plt.subplots(figsize=(6.6, 4.9))

for name, spec in constructions.items():
    F = fidelity(T1_ns_list, spec["N_pulses"])
    infid = 1 - F
    ax.plot(T1_us_list, infid, color=spec["color"], lw=1.8,
            marker=spec["marker"], markevery=20, ms=6,
            label=fr"{name} ($N={spec['N_pulses']}$)")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$T_1$ ($\mu$s)")
ax.set_ylabel(r"Projected Toffoli infidelity, $1-F$")
ax.set_title("Projected Toffoli-gate infidelity vs. $T_1$\n"
             r"(common per-pulse duration $t_{\rm pulse}=132.75$ ns)", fontsize=12)
ax.legend(loc="lower left", fontsize=8.7, frameon=True)
ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.6)

fig.tight_layout()
fig.savefig("toffoli_noise_comparison_v2.png", dpi=220)
print("Saved toffoli_noise_comparison_v2.png")

# ----------------------------------------------------------------------
# 4. Print representative values for the write-up
# ----------------------------------------------------------------------
for T1_us in [10, 100, 1000]:
    T1_ns = T1_us * 1e3
    print(f"\nT1 = {T1_us} us:")
    for name, spec in constructions.items():
        F = fidelity(T1_ns, spec["N_pulses"])
        print(f"  {name:24s} N={spec['N_pulses']:2d}  F={F:.4f}")
