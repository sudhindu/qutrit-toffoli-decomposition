# Qutrit-Assisted Toffoli / C^k-NOT Decomposition — Simulation Code

This repository contains the complete simulation code, generated
data, and figures supporting the numerical validation sections of the
paper *"[paper title]"* (Sudhindu Bikash Mandal, CQuERE). It covers
the full chain of simulations: single-qutrit gate validation, the
corresponding noise analysis, two-qutrit controlled-gate realization
via cross-resonance driving, and a projected full-circuit noise
comparison against three prior qutrit-based constructions.

**Every script has its own dedicated documentation page** in
[`docs/`](docs/), covering the physical theory, the algorithm, and —
critically — an explicit account of which parameters were
*numerically derived or calibrated* versus *chosen as representative
operating points*, so that no parameter in this codebase is presented
as more justified than it actually is.

---

## Repository structure

```
qutrit-toffoli-decomposition/
├── README.md                              <- you are here
├── requirements.txt
├── docs/
│   ├── 01_single_qutrit_gates.md          <- theory + algorithm for X01/X12
│   ├── 02_single_qutrit_noise.md          <- theory + algorithm for the noise sweep
│   ├── 03_two_qutrit_controlled_gates.md  <- theory + algorithm for C1X12/C2X01/C2X12
│   └── 04_toffoli_noise_projection.md     <- theory + algorithm for the Toffoli-level comparison
├── src/
│   ├── single_qutrit_gates/
│   │   ├── transmon_x01_x12_sim.py
│   │   └── transmon_x01_x12_noise_sim.py
│   ├── two_qutrit_controlled_gates/
│   │   └── transmon_C1X12_C2X01_C2X12_CR_sim_v3.py
│   └── toffoli_noise_projection/
│       └── toffoli_noise_comparison_v2.py
└── results/
    ├── single_qutrit_gates/               <- figures, CSVs, and text tables already generated
    ├── two_qutrit_controlled_gates/
    └── toffoli_noise_projection/
```

---

## What each module validates, in one line

| # | Module | Validates | Doc page |
|---|---|---|---|
| 1 | `single_qutrit_gates` (ideal) | $X_{01}$, $X_{12}$ single-qutrit gates, via full lab-frame (no-RWA) simulation | [`docs/01`](docs/01_single_qutrit_gates.md) |
| 2 | `single_qutrit_gates` (noise) | Same gates under Lindblad relaxation/dephasing vs. $T_1$ | [`docs/02`](docs/02_single_qutrit_noise.md) |
| 3 | `two_qutrit_controlled_gates` | $C^1X_{12}$, $C^2X_{01}$, $C^2X_{12}$ via cross-resonance driving between two coupled transmons | [`docs/03`](docs/03_two_qutrit_controlled_gates.md) |
| 4 | `toffoli_noise_projection` | Projected full-Toffoli fidelity vs. $T_1$, comparing this work's pulse count against three prior qutrit-based constructions | [`docs/04`](docs/04_toffoli_noise_projection.md) |

---

## Quick start

```bash
git clone <this-repo-url>
cd qutrit-toffoli-decomposition
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -r requirements.txt

# Run any module, e.g.:
python src/single_qutrit_gates/transmon_x01_x12_sim.py
```

Each script writes its figures/CSVs to its own working directory by
default — see the "How to run" section at the bottom of each doc page
in `docs/` for expected runtime and output filenames. The `results/`
folder already contains the exact output artifacts referenced in the
paper, so re-running the scripts is only necessary to reproduce or
extend them, not to view the reported results.

---

## Key methodological points (see `docs/` for full derivations)

- **No rotating-wave pre-approximation.** All Hamiltonians are
  integrated in the lab frame, carrier included — every reported
  leakage/infidelity number is a genuine simulated result, not an
  artifact of an analytic approximation baked into the model.
- **Every calibrated parameter is derived, not guessed.** Single-qutrit
  pulse amplitudes are solved exactly from the analytic $\pi$-pulse
  area condition; two-qutrit controlled-gate amplitudes *and*
  durations are found via an explicit numerical joint-optimization
  scan (no closed-form solution exists for that coupled system) — see
  [`docs/03`](docs/03_two_qutrit_controlled_gates.md) §4 for the exact
  procedure.
- **A documented failure mode is kept in the repository, not erased.**
  An earlier, incorrect choice of target transmon frequency produced
  misleading simulated results due to an accidental spectral
  near-degeneracy; the diagnostic that caught this (and must be re-run
  whenever device parameters change) is built into the two-qutrit
  script and explained in [`docs/03`](docs/03_two_qutrit_controlled_gates.md) §3.
- **Weak/imperfect results are reported as-is.** The $C^1X_{12}$
  gate's weakest control-selectivity branch (69.3%) is reported
  directly in the output table, not adjusted or omitted — see
  [`docs/03`](docs/03_two_qutrit_controlled_gates.md) §6.
- **The Toffoli-level comparison (`toffoli_noise_comparison_v2.py`)
  contains no quantum simulation at all** — it is an explicitly scoped,
  literature-grounded projection, not a re-simulation of competing
  hardware. See [`docs/04`](docs/04_toffoli_noise_projection.md) §1 for
  exactly what it does and does not claim.

---

## Requirements

See [`requirements.txt`](requirements.txt). Tested with Python 3.11–3.13
and QuTiP 5.x.

## Citation

If you use this code, please cite:

```
[BibTeX entry for the paper — add once available]
```

## License

[Add your chosen license here, e.g. MIT — see `LICENSE`]
