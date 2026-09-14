# `transmon_x01_x12_sim.py` — Single-Qutrit Gate Validation

**Location:** `src/single_qutrit_gates/transmon_x01_x12_sim.py`
**Outputs:** `results/single_qutrit_gates/X01_gate_simulation.png`,
`X12_gate_simulation.png`, `transfer_tables.txt`

---

## 1. What this script validates

The paper's Toffoli/$C^k$-NOT construction is built out of two
single-qutrit primitives, $X_{01}$ and $X_{12}$, each defined to act
as a NOT gate on one two-level subspace of a three-level qutrit while
leaving the third level completely untouched:

$$X_{01} = |1\rangle\langle 0| + |0\rangle\langle 1| + |2\rangle\langle 2|, \qquad
  X_{12} = |0\rangle\langle 0| + |2\rangle\langle 1| + |1\rangle\langle 2|.$$

This script asks a specific, falsifiable question: **if a real
microwave pulse is applied to a real (simulated) transmon Hamiltonian,
does the resulting physical evolution actually reproduce these ideal
matrices, to what precision, and is that precision a genuine physical
result rather than an assumption baked into the model?**

---

## 2. Physical model (theory)

### 2.1 Bare transmon Hamiltonian

The qutrit is modeled as an anharmonic oscillator with three levels
$|0\rangle,|1\rangle,|2\rangle$. Setting $\hbar=1$ and the ground-state
energy to zero:

$$H_0 = \omega_{01}\,|1\rangle\langle1| + (2\omega_{01}-\alpha)\,|2\rangle\langle2|$$

where $\omega_{01}$ is the $0\!\leftrightarrow\!1$ transition
frequency and $\alpha$ is the anharmonicity. Note the $|2\rangle$
energy is $2\omega_{01}-\alpha$, **not** $2\omega_{01}$: a perfectly
harmonic ladder would have $\alpha=0$; the anharmonicity is precisely
what makes the third level distinguishable and separately addressable
from the second.

### 2.2 Drive Hamiltonian

A microwave drive couples to the qutrit's charge operator. The
$\sqrt{2}$ factor on the $1\!\leftrightarrow\!2$ term is not a free
parameter — it is the standard harmonic-oscillator-ladder scaling of
the dipole matrix element for the $n\!\leftrightarrow\!(n+1)$
transition ($\langle n+1|\hat a^\dagger|n\rangle \propto \sqrt{n+1}$):

$$H_d(t) = \Omega(t)\cos(\omega_d t)\Big(|0\rangle\langle1|+|1\rangle\langle0|
           + \sqrt2\,(|1\rangle\langle2|+|2\rangle\langle1|)\Big)$$

### 2.3 Why no rotating-wave approximation (RWA)

The script integrates $H(t)=H_0+H_d(t)$ **exactly as written**, with
the fast carrier $\cos(\omega_d t)$ kept explicit — it does not
transform into a rotating frame or discard counter-rotating terms
before simulating. This is a deliberate methodological choice: the
paper's analytic derivation of gate isolation (Sec. II of the paper)
*does* use the RWA to obtain a clean effective two-level Hamiltonian.
Simulating the **un-approximated** lab-frame equation is what lets
this script serve as an independent numerical check on that
analytic approximation, rather than a restatement of it — if the RWA
had discarded something non-negligible, this simulation is exactly
where that would show up as unexpected leakage.

---

## 3. Algorithm

```
INPUT:  f01, alpha           (transmon frequencies, GHz)
        T_PULSE, N_STEPS     (pulse duration and time resolution)

1.  Build H0 from (f01, alpha)                        [Sec 2.1]
2.  Build the fixed drive operator Hd_op = b01 + sqrt(2) b12
3.  FOR each gate in {X01, X12}:
        a. Choose carrier frequency wd
             X01: wd = w01           (resonant with 0<->1)
             X12: wd = w12 = w01-alpha (resonant with 1<->2)
        b. Build a Gaussian pulse envelope Omega(t), centered in the
           pulse window, width = T_PULSE/6
        c. SOLVE for the envelope amplitude such that the pulse AREA
           exactly satisfies the analytic pi-pulse condition
             X01:  integral( Omega(t) dt )         = pi
             X12:  integral( sqrt(2) Omega(t) dt )  = pi
        d. Assemble H(t) = H0 + Omega(t) cos(wd t) Hd_op
        e. FOR each initial state psi0 in {|0>, |1>, |2>}:
               integrate the time-dependent Schrodinger equation
               (qutip.sesolve) over [0, T_PULSE]
               record P0(t), P1(t), P2(t) at every time step
        f. Assemble the 3x3 final-time population-transfer matrix
4.  Plot pulse envelope + population dynamics (4-panel figure)
5.  Save the two 3x3 transfer-table matrices to transfer_tables.txt
```

---

## 4. Parameter justification — what is calibrated vs. what is representative

This is the question the paper's reviewers (and you) should be able
to answer precisely for every number in the script — the table below
states, for each parameter, whether it was **derived** from a stated
physical condition, or **chosen as a representative value** (and if
so, why that choice is defensible).

| Parameter | Value | Status | Justification |
|---|---|---|---|
| $f_{01}$ | 5.000 GHz | Representative | A standard, realistic transmon 0–1 frequency; not fit to any specific device, chosen to be in the range reported across the comparison literature (Gokhale, Galda, Nikolaeva). |
| $\alpha$ | 300 MHz | Representative | Typical fixed-frequency transmon anharmonicity magnitude (order 200–350 MHz is standard in the literature). |
| $\Omega_{\max}$ (peak Rabi rate) | ~30.00 MHz ($X_{01}$), ~21.21 MHz ($X_{12}$) | **Derived, not assumed** | These are *not* independently chosen — they are the unique amplitudes that make the discretized pulse area exactly equal to the analytic $\pi$ (or $\pi/\sqrt2$) condition, for the *given* envelope shape and duration (Sec. 5 below). Changing the shape or duration changes this number; it is a solved quantity, not a free input. The exact ratio between the two ($1/\sqrt2$) is a direct, checkable consequence of the $\sqrt2$ dipole scaling in §2.2 — this ratio coming out correctly in the simulation is itself a validation that the code is self-consistent with the stated theory. |
| Pulse shape | Truncated Gaussian | Chosen (standard practice) | A smooth, rather than abrupt, envelope is standard experimental practice — it avoids injecting spectral weight into unwanted frequencies (a rectangular pulse's sharp edges correspond to a broad frequency content). This is a *shape* choice, not a fitted physical parameter. |
| $T_{\text{PULSE}}$ | 40 ns | Representative, **not optimized** | This is the one parameter in this script that is genuinely a chosen operating point rather than a derived or scanned quantity: 40 ns is a realistic, commonly reported single-qutrit gate duration for transmons, but no scan over duration was performed here to find an "optimal" value. (Contrast this with the two-qutrit controlled-gate script, `transmon_C1X12_C2X01_C2X12_CR_sim_v3.py`, where both amplitude *and* duration genuinely are the result of a numerical scan — see that script's doc page.) The honest characterization of this script's parameters is: two physical constants chosen to be realistic (not fit), one shape chosen for standard practice, one duration chosen as a representative operating point, and the amplitude solved exactly from the other four. |
| $N_{\text{STEPS}}$ | 40,000 (dt = 1 ps) | Derived from a resolution requirement | Not arbitrary: the carrier oscillates at $f_{01}=5$ GHz, i.e. one cycle every 200 ps. 1 ps steps give ~200 samples per carrier cycle, comfortably above the minimum needed to numerically resolve the fastest oscillation in the problem without aliasing. |

---

## 5. How the pulse amplitude is actually solved (not assumed)

```python
def build_envelope(area_target, coupling_factor=1.0):
    shape = gaussian_shape(tlist, t0, sigma)              # unit-height shape
    raw_area = np.trapezoid(shape, tlist)                 # numeric area of that shape
    amplitude = area_target / (coupling_factor * raw_area)
    return amplitude * shape
```

This is solving one linear equation:
`coupling_factor * amplitude * raw_area = area_target`. Given the
*shape* (fixed once $T_{\text{PULSE}}$ and $\sigma$ are chosen) and
the *target area* (fixed by the analytic $\pi$-pulse condition), the
amplitude is the **unique** number satisfying both — it is a solved
output of the calibration condition, not a separately guessed input.

---

## 6. Validation checks built into the output

- **Population conservation:** $P_0(t)+P_1(t)+P_2(t)=1$ at every
  instant (a basic sanity check available by inspecting the plotted
  curves).
- **Spectator invariance:** the level the gate is *not* supposed to
  touch stays at population 1.000 throughout — this is the main
  physical claim being tested, and it is a genuine output of the
  unapproximated simulation, not an assumption.
- **Forward/reverse symmetry:** $X_{01}$ maps $|0\rangle\!\to\!|1\rangle$
  and $|1\rangle\!\to\!|0\rangle$ with matching fidelity (99.8% both
  directions in the reported run) — consistent with $X_{01}$ being a
  true involution ($X_{01}^2=I$) at the level of a clean two-level
  Rabi problem.

---

## 7. How to run

```bash
python transmon_x01_x12_sim.py
```

Requires: `numpy`, `qutip`, `matplotlib`. Runtime: a few seconds on a
standard laptop (6 `sesolve` calls total, each over a 40 ns / 40,000-step
window).
