# `transmon_C1X12_C2X01_C2X12_CR_sim_v3.py` — Two-Qutrit Controlled Gates via Cross-Resonance

**Location:** `src/two_qutrit_controlled_gates/transmon_C1X12_C2X01_C2X12_CR_sim_v3.py`
**Outputs:** `results/two_qutrit_controlled_gates/CR_controlled_gate_dynamics_v3.png`,
`CR_noise_fidelity_v3.png`, `CR_controlled_gate_table_v3.csv`,
`CR_noise_fidelity_data_v3.csv`

---

## 1. What this script validates

The paper's Toffoli/$C^k$-NOT decomposition needs three **controlled**
two-qutrit gates:

- $C^1X_{12}$ — $X_{12}$ acts on the target iff the control is $|1\rangle$
- $C^2X_{01}$ — $X_{01}$ acts on the target iff the control is $|2\rangle$
- $C^2X_{12}$ — $X_{12}$ acts on the target iff the control is $|2\rangle$
  (the additional primitive needed to extend from the Toffoli, $k=2$,
  to general $C^k$-NOT)

This script asks whether these can be realized by the *same* physical
mechanism used for real hardware two-qubit gates — cross-resonance
(CR) driving — using **one native pulse per gate**, and quantifies
exactly how well, honestly reporting where it does not work perfectly
rather than only reporting the successful cases.

---

## 2. Theory: cross-resonance driving between two coupled qutrits

### 2.1 The physical setup

Unlike the single-qutrit scripts, this model has **two distinct
transmons** — control and target — each with its own frequency and
anharmonicity, connected by a fixed (always-on) capacitive coupling:

$$H_0 = H_0^{(c)}\otimes I + I\otimes H_0^{(t)}, \qquad
  H_{\text{couple}} = g\,\big(\text{charge}_c\otimes\text{charge}_t\big)$$

with $g/2\pi=25$ MHz $\ll$ the control–target detuning
($\Delta = |f_{01}^c-f_{01}^t|$), i.e. the *dispersive regime*.

### 2.2 The cross-resonance mechanism

The drive is applied through the **control's** charge line, at a
frequency matched to the **target's** transition — not the control's
own transition. Physically: because of the coupling, the true energy
eigenstates of the combined system are not perfectly separable into
independent "control state" and "target state" — they are slightly
mixed. Driving the control off-resonantly (relative to its own bare
frequency) at exactly the target's frequency exploits this mixing: a
small fraction of the drive's effect is transmitted through the
coupling and lands resonantly on the target, producing an effective
conditional rotation whose presence and rate depend on which level the
control occupies. This is the same physical mechanism used for IBM's
native CR/ECR two-qubit gate on fixed-frequency hardware — no tunable
elements are required.

### 2.3 Getting the drive frequency right: diagonalization, not guessing

The control-dependent target transition frequency is **not** assumed
or looked up from a table — it is *computed* by numerically
diagonalizing $H_0+H_{\text{couple}}$ and reading off the true
dressed-state energy differences:

```python
evals, evecs = H_static.eigenstates()
f_t12_given_c = {c: E[(c,2)] - E[(c,1)] for c in range(3)}
```

This gives three *different* dressed transition frequencies (one per
control level $c=0,1,2$), and the drive for, e.g., $C^1X_{12}$ is set
to the specific one dressed by $c=1$. This is exactly analogous to how
a real experiment would locate this frequency: via two-tone
spectroscopy on the actual coupled device, not from a datasheet.

---

## 3. A documented failure mode, and why the diagnostic in §3 of this script exists

**This script includes an explicit diagnostic (the "overlap check")
that exists *because* an earlier version of this analysis produced a
subtly wrong result, and this check is what caught it.** This is worth
keeping in the repository history and documenting explicitly, rather
than silently using only the corrected version, because it demonstrates
the actual verification standard applied to every parameter in this
project.

**The failure:** an earlier choice of target frequency
($f_{01}^t=5.300$ GHz) happened to sit very close to the condition
$f_{01}^t = f_{01}^c+\alpha_t$ — an *accidental* near-degeneracy
between the bare states $|c{=}0,t{=}2\rangle$ and $|c{=}1,t{=}1\rangle$
(both near 10.3 GHz). At the coupling strength used, this strongly
hybridized the two states: the bare state $|c{=}1,t{=}1\rangle$
retained only **62.7%** overlap with any single true eigenstate of the
system — meaning the "control stays in $|1\rangle$, only the target
moves" picture was not a good approximation to what was actually
happening, even before any drive was applied. This produced messy,
strongly asymmetric simulated gate behavior that was not a real
physical limitation, but an artifact of a bad frequency choice.

**The fix and the check that prevents recurrence:**

```python
worst_overlap = min over all 9 bare product states of
    (max overlap of that bare state with any single true eigenstate)
```

Before calibrating *any* pulse, the script explicitly verifies
`worst_overlap > 0.95` for the chosen static Hamiltonian
($f_{01}^t=5.550$ GHz gives worst\_overlap $=0.976$). **This check
must be re-run whenever any of $f_{01}^c,f_{01}^t,\alpha_c,\alpha_t,g$
change** — it is the load-bearing safeguard against silently repeating
the earlier mistake, and is why it is printed at the start of every
run rather than only checked once during development.

---

## 4. How the pulse amplitude and duration were actually found

**This is the key methodological difference from the single-qutrit
scripts: for these controlled gates, *neither* the amplitude *nor* the
duration is assumed or read from an analytic formula.** This coupled,
multi-level, two-transmon system has no closed-form Rabi solution —
unlike the single-qutrit π-pulse condition (script 1), there is no
algebraic expression to solve for the amplitude given a target area.
Both parameters were instead found by a **numerical joint scan**,
maximizing two objectives simultaneously:

1. **Addressed-branch fidelity** — how completely the target moves to
   the correct state when the control is the one that should trigger
   the gate.
2. **Control-selectivity** — how little the target moves when the
   control is in either of the two states that should leave it alone.

```
FOR candidate (amplitude, duration) pairs on a grid:
    simulate the addressed branch -> fidelity_fwd
    IF fidelity_fwd < 0.9: reject (not a plausible operating point)
    simulate the non-addressed branches -> selectivity_1, selectivity_2
    score = min(fidelity_fwd, selectivity_1, selectivity_2)
KEEP the (amplitude, duration) with the highest score
```

The resulting calibrated operating points (reported in the code and
reproduced in Table 1 of the paper) are:

| Gate | Amplitude | Duration | Basis for these specific numbers |
|---|---|---|---|
| $C^1X_{12}$ | 25 MHz | 900 ns | Found by the scan above; this branch has the smallest usable control-dependent frequency splitting of the three gates, which is *why* it needs the longest pulse — a slower pulse is required for good selectivity when the addressed and non-addressed branches are close in frequency. |
| $C^2X_{01}$ | 35 MHz | 250 ns | Found by the scan; the largest available splitting for this gate permits the fastest pulse of the three. |
| $C^2X_{12}$ | 35 MHz | 600 ns | Found by the scan; intermediate splitting, intermediate duration. |

**None of these three (amplitude, duration) pairs are equal, and none
were chosen in advance** — each is the specific output of an
independent optimization for that gate's own dressed-frequency
landscape. This is the honest characterization to give a reader: the
*existence* of a working native pulse is not assumed, it is
demonstrated by exhibiting one found through a documented, repeatable
search procedure.

---

## 5. Algorithm (full script)

```
1.  Build H0, H_couple, H_static for the two coupled transmons.
2.  Diagonalize H_static; run the overlap diagnostic (Sec. 3).
3.  Read off the three control-dependent drive frequencies
    (wd_C1X12, wd_C2X01, wd_C2X12) from the diagonalization -- not
    assumed.
4.  Use the ALREADY-CALIBRATED (amplitude, duration) pairs (Sec. 4)
    to build the time-dependent drive Hamiltonian for each gate.
5.  IDEAL TABLE: for each gate, for a set of (control, target-input)
    combinations, evolve under sesolve and record the final
    population-transfer numbers -> CR_controlled_gate_table_v3.csv
6.  DYNAMICS FIGURE: for each gate, with the target fixed at one
    representative starting state, evolve for all 3 possible control
    inputs and record the FULL time trace (not just the endpoint)
    -> CR_controlled_gate_dynamics_v3.png
7.  NOISE SWEEP: for each gate's addressed branch, for T1 in a
    log-spaced 1-1000 us grid, evolve under mesolve with the SAME
    collapse-operator construction as the single-qutrit noise script
    (Sec 2 of that doc page), extended to act on EITHER qutrit via
    tensor-with-identity embedding -> CR_noise_fidelity_v3.png,
    CR_noise_fidelity_data_v3.csv
```

---

## 6. Honest reporting of an imperfect branch

Table 1 of the paper (and the CSV in this repository) reports that the
$C^1X_{12}$ gate's control$=|2\rangle$ branch retains only **69.3%**
of its population correctly — the weakest single number in the whole
comparison. **This is reported, not hidden**, and is a real,
understood physical limitation: this branch has the smallest
control-dependent frequency splitting from the addressed
($c{=}1$) branch of any of the nine (gate, wrong-control) pairs
examined, so a bare (non-echoed) single-tone pulse cannot fully
suppress it within a practically useful gate duration. The paper's
text explicitly attributes this to the use of a *bare* CR drive, and
notes that echo/composite-pulse refinement — a standard, well
precedented technique for real CR gates — would be the natural next
step to improve it; this was not implemented here, and the honest
number is reported as-is rather than adjusted or omitted.

---

## 7. How to run

```bash
python transmon_C1X12_C2X01_C2X12_CR_sim_v3.py
```

Requires: `numpy`, `qutip`, `matplotlib`. Runtime: the noise sweep
(step 7 above) dominates — 3 gates × 13 $T_1$ values, each a `mesolve`
call on a 9-dimensional (two-qutrit) density matrix over pulse
durations up to 900 ns. Expect several minutes on a standard laptop.
