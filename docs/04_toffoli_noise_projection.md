# `toffoli_noise_comparison_v2.py` — Projected Full-Toffoli Fidelity Comparison

**Location:** `src/toffoli_noise_projection/toffoli_noise_comparison_v2.py`
**Output:** `results/toffoli_noise_projection/toffoli_noise_comparison_v2.png`

---

## 1. What this script is — and, importantly, what it is *not*

This script is fundamentally different in kind from the other three:
it does **not** simulate any Hamiltonian. It contains no `qutip` call
at all. This is a deliberate scope decision, stated explicitly here so
it is never mistaken for something it is not: fully re-simulating the
*other three* papers' Toffoli constructions would require independently
modeling each paper's own native two-qutrit gate (Nikolaeva's native
iSWAP⁰², Galda's echoed cross-resonance CNOT reused via an incidental
rotation-halving trick, Gokhale's cyclic-gate decomposition), each of
which is its own separate hardware/theory problem outside the scope of
validating *this paper's* construction. Instead, this script builds a
**transparent, literature-grounded, gate-count-based projection** —
standard practice for this kind of cross-construction comparison —
that isolates the effect of one specific, already-established
quantity (total physical pulse count per construction) on projected
fidelity, holding everything else fixed.

---

## 2. Theory: the idle-limited average gate infidelity approximation

For a single pulse of duration $t$ acting on a system with relaxation
time $T_1$, the standard small-$t/T_1$ approximation for average gate
infidelity under pure amplitude damping is:

$$1-F_{\text{pulse}} \approx \frac{2}{3}\cdot\frac{t}{T_1}$$

This is a widely used, standard formula in the gate-fidelity/randomized-benchmarking
literature — not a *bespoke* choice for this comparison. For a circuit
of $N$ **independent** sequential pulses, each contributing this same
per-pulse infidelity, the *combined* fidelity of the full circuit is
the product:

$$F_{\text{total}}(T_1) = \Big(1-\tfrac23\tfrac{t_{\text{pulse}}}{T_1}\Big)^{N}$$

This is the single formula the whole script evaluates, once per
construction, varying only $N$ (pulse count) and — held fixed across
all constructions — $t_{\text{pulse}}$ (per-pulse duration).

---

## 3. Where the common per-pulse duration comes from — grounded in a real number, not guessed

$$t_{\text{pulse}} = \frac{1.593\ \mu\text{s}}{12} = 132.75\ \text{ns}$$

The numerator, $1.593\ \mu$s, is **Galda et al.'s own reported total
experimental gate time** for their ternary CCNOT decomposition on
`ibmq_jakarta` — a real, independently citable measured number, not an
assumption invented for this comparison.

The denominator, **12**, is the *audited* total physical pulse count
for that same construction (paper §III, `docs/toffoli_subsection`
derivation): two $C^1X_{12}$-equivalent fold/unfold blocks at 5 pulses
each ($X_-$: 2 pulses, CNOT: 1 pulse, $X_+$: 2 pulses), plus their
middle $|2\rangle$-controlled target flip at 2 pulses (their own stated
realization: reusing a single control$=|1\rangle$-calibrated CNOT pulse
twice) — $5+5+2=12$. **This replaces an earlier, incorrect version of
this calculation that divided by 4** (the entangling-pulse-only count
originally reported in the comparison table), which implied an
inflated $t_{\text{pulse}}=398.25$ ns; using the correct 12-pulse
denominator gives the corrected, smaller $132.75$ ns figure used here.
This correction is documented so that anyone re-deriving this number
from the cited paper arrives at the same 12, not the earlier 4.

**This one real number is then applied uniformly to every
construction in the comparison** — including this paper's own
3-pulse Toffoli — deliberately, so that the comparison isolates the
effect of *pulse count* specifically, rather than each paper's own
device-specific pulse speed (which was not independently re-derived
for the other constructions, and is explicitly out of scope per §1).

---

## 4. Algorithm

```
INPUT:  pulse counts N, one per construction (from the paper's Table 1):
            This work                 N = 3
            This work + Galda's CX    N = 8   (hybrid; see paper)
            Gokhale et al.            N = 6
            Galda et al. (audited)    N = 12
            Nikolaeva et al.          N = 4
        t_pulse = 132.75 ns            (Sec. 3, common to all)

FOR T1 in 200 log-spaced points from 1 us to 1000 us:
    FOR each construction (name, N):
        F(T1) = (1 - (2/3)*(t_pulse/T1))^N
        infidelity(T1) = 1 - F(T1)

Plot infidelity vs T1 (log-log) for all constructions on one figure.
Print representative F(T1) values at T1 = 10, 100, 1000 us.
```

---

## 5. Parameter justification summary

| Parameter | Value | Status |
|---|---|---|
| $t_{\text{pulse}}$ | 132.75 ns | **Derived** from a real, cited experimental number (Galda et al.'s reported total gate time) divided by an independently *audited* pulse count — not fit or guessed, and explicitly corrected once already (Sec. 3). |
| Pulse counts $N$ (all 5 constructions) | 3, 4, 6, 8, 12 | **Derived**, each via an explicit pulse-by-pulse decomposition documented in the paper's Toffoli subsection — none are taken uncritically from the other papers' own summary tables (the Galda number specifically was found to have been mislabeled in an earlier draft of this comparison, and was corrected here). |
| $\frac{2}{3}(t/T_1)$ formula | — | **Standard**, cited approximation from the gate-fidelity literature, not custom-derived for this comparison. |
| $T_1$ range | 1–1000 μs | Chosen to span the same realistic device range used throughout this repository's other noise sweeps, for visual/narrative consistency across all figures in the paper. |

---

## 6. What this script deliberately does not claim

- It does **not** claim that the other three constructions would
  achieve exactly these numbers on real hardware — only that *if* each
  physical pulse in every construction took the same characteristic
  duration (grounded in one real reported number), the resulting
  fidelity gap would be driven entirely by pulse count, which is
  precisely the quantity this paper's construction reduces.
- It does **not** re-derive Gokhale's or Nikolaeva's pulse counts from
  first principles the way it does for Galda's (Sec. 3) — those two
  entries are carried over from the paper's existing, separately
  justified Table 1 values. A reader auditing this script should treat
  the Galda number as independently re-derived here, and the other two
  as sourced from the paper's main comparison table.

---

## 7. How to run

```bash
python toffoli_noise_comparison_v2.py
```

Requires: `numpy`, `matplotlib` only (no `qutip` — see Sec. 1).
Runtime: instantaneous (a closed-form formula evaluated on a fixed
grid, no numerical integration).
