# `transmon_x01_x12_noise_sim.py` — Single-Qutrit Gate Under Realistic Decoherence

**Location:** `src/single_qutrit_gates/transmon_x01_x12_noise_sim.py`
**Outputs:** `results/single_qutrit_gates/X01_X12_noise_fidelity.png`,
`noise_fidelity_data.csv`
**Depends on:** the same physical model as
[`01_single_qutrit_gates.md`](01_single_qutrit_gates.md) — read that
page first; this page covers only what is *new* here (the noise
model).

---

## 1. What this script adds

The previous script simulates a **closed** quantum system — no
interaction with the environment. Real superconducting qubits are
never perfectly isolated. This script asks: **how does gate fidelity
degrade as a function of realistic energy-relaxation and dephasing
timescales, and does the degradation vanish correctly in the limit of
a perfectly isolated qubit ($T_1\to\infty$)?** The second half of that
question is a built-in, non-trivial self-consistency check (Sec. 4).

---

## 2. Theory: the Lindblad master equation

A closed system evolves under the Schrödinger equation for a state
vector $|\psi\rangle$. An **open** system — one that exchanges energy
and information with an uncontrolled environment — must instead be
described by a density matrix $\rho$ (which can represent statistical
mixtures, not just pure superpositions), evolving under the Lindblad
master equation:

$$\dot\rho = -i[H(t),\rho] + \sum_k \Big(L_k \rho L_k^\dagger
              - \tfrac12\{L_k^\dagger L_k,\rho\}\Big)$$

Each $L_k$ (a "collapse operator" or "jump operator") represents one
physical decoherence channel. This script defines three:

| Collapse operator | Physical process | Rate |
|---|---|---|
| $L_1=\sqrt{\Gamma_1}\,|0\rangle\langle1|$ | Spontaneous relaxation $|1\rangle\to|0\rangle$ | $\Gamma_1=1/T_1$ |
| $L_2=\sqrt{2\Gamma_1}\,|1\rangle\langle2|$ | Spontaneous relaxation $|2\rangle\to|1\rangle$ | $2\Gamma_1$ |
| $L_3=\sqrt{2\Gamma_1}\,\mathrm{diag}(0,1,2)$ | Pure dephasing (no energy change, phase randomization) | $2\Gamma_1$, with $T_\phi\equiv T_1$ |

**The factor of 2 on $L_2$ is not an independent assumption** — it is
required for internal consistency with the coherent drive operator
already defined in the closed-system model
($H_d \propto |1\rangle\langle2|\sqrt2 + \text{h.c.}$): the same
$\sqrt{n+1}$ ladder-operator scaling that enhances the *drive*
coupling to the $1\!\leftrightarrow\!2$ transition by $\sqrt2$ also
enhances the *spontaneous emission* rate on that same transition by
$(\sqrt2)^2=2$ — both effects trace back to the same physical dipole
matrix element. Choosing an independent, unrelated relaxation rate for
$|2\rangle\to|1\rangle$ would have been the actual unjustified
assumption; this factor of 2 is forced by the model already committed
to elsewhere.

**Setting $T_\phi=T_1$** (pure dephasing rate equal to relaxation
rate) is an explicit simplification, stated as such: it collapses two
independent noise timescales into one swept parameter ($T_1$) in the
absence of an independently measured $T_2$. This is a standard
default in the noise-modeling literature when only one coherence
timescale is available or assumed.

---

## 3. The fidelity metric: why "truth-table fidelity," not QuTiP's built-in process fidelity

$$F_{TT} = \frac{1}{3}\sum_{i=0}^{2} \langle\text{target}_i|\,\rho_{\text{out}}(\psi_0=|i\rangle)\,|\text{target}_i\rangle$$

— the population found in the *correct* output state, averaged over
the three relevant computational-basis inputs. This choice was made
deliberately, not as a default: during development, QuTiP's
higher-level Choi-matrix/process-fidelity routines were found to
disagree with the direct, hand-computed truth-table fidelity for this
specific time-dependent, non-rotating-frame Hamiltonian (traced to the
propagator construction, not the physics). $F_{TT}$ was adopted
because every quantity feeding into it can be independently checked by
eye against the same population curves produced in the first script —
it is the more *auditable* metric, which matters more here than using
a nominally more sophisticated one that could not be independently
verified.

---

## 4. Algorithm

```
INPUT: T1_us_list = 16 points, log-spaced from 1 us to 1000 us

1.  Build H(t) for X01 and X12 exactly as in the noiseless script.
2.  self-consistency check:
        F_check = truth_table_fidelity(H, targets, c_ops=[])
        (i.e. call the NOISY solver with an EMPTY collapse-operator
         list, and confirm it reproduces the independently-computed
         noiseless result from script 1: 99.88% / 99.99%)
3.  FOR each T1 in T1_us_list:
        build collapse_ops(T1)
        FOR each gate in {X01, X12}:
            FOR each of the 3 basis-state inputs:
                rho0 = |psi0><psi0|
                evolve rho0 under mesolve(H, rho0, c_ops, [0,T_PULSE])
                record population in the correct target state
            F_TT(T1) = average of the 3 populations
4.  Plot infidelity (1 - F_TT) vs T1 on a log-log scale, for both
    gates, with the T1->infinity noiseless floor drawn as a
    horizontal reference line for each.
5.  Save (T1, F_TT_X01, F_TT_X12, infidelities) to CSV.
```

---

## 5. Parameter justification

| Parameter | Value | Status | Justification |
|---|---|---|---|
| $T_1$ sweep range | 1–1000 μs, 16 log-spaced points | Chosen to span the literature | Covers early-generation transmons (~1–10 μs) through state-of-the-art devices (100s of μs–1 ms) reported in the comparison papers, so the plot is informative across the full range a reader might care about — not fit to any single target device. |
| $T_\phi = T_1$ | — | Explicit simplification | See §2 above; stated, not hidden. |
| All physical Hamiltonian parameters ($f_{01},\alpha,\Omega_{\max}$, pulse shape/duration) | — | **Identical to script 1** | Deliberately unchanged, so that any difference between the two scripts' results is attributable *only* to the newly introduced noise, not to a re-tuned Hamiltonian. |

---

## 6. The built-in self-consistency check, explained

```python
f_check_X01 = truth_table_fidelity(H_X01, targets_X01, [])
```

Passing an *empty* collapse-operator list to the noisy solver
(`mesolve`) reduces the Lindblad equation exactly to the Schrödinger
equation — this is a mathematical identity, not an approximation. So
this single line re-derives the noiseless result using the entire
noisy-solver code path, and the script explicitly checks it against
the *independently* computed result from script 1
(99.88% / 99.99%). If these two numbers had disagreed, it would
indicate a bug in the noise machinery, independent of whether the
underlying physics was correct — this is why the check exists, and
why passing it is reported as meaningful evidence, not a formality.

The same check reappears in the **T1 → ∞ limit** of the actual sweep:
at $T_1=1000\ \mu$s, $F_{TT}$ should closely approach (not necessarily
exactly reach, since 1000 μs is large but finite) the same 99.88%/99.99%
floor — and the plotted dotted reference lines make this visually
checkable in the figure itself.

---

## 7. How to run

```bash
python transmon_x01_x12_noise_sim.py
```

Requires: `numpy`, `qutip`, `matplotlib`. Runtime: this script calls
`mesolve` $3\text{ inputs}\times16\ T_1\text{ values}\times2\text{ gates}=96$
times — expect roughly 1–2 minutes on a standard laptop, noticeably
longer than script 1.
