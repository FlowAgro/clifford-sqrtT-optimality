# Certified optimality for Clifford+√T gate synthesis

Independent, exhaustive-search measurements of how close published Clifford+√T single-qubit circuits are to the
true minimum cost, and of optimal randomized (mixed) synthesis in this gate set, with every central result
re-verified through an independent code path.

Author: **Cristian Vitorino Rodrigues da Cunha**, independent researcher, Uberlândia, Brazil
(ORCID [0009-0009-4465-6005](https://orcid.org/0009-0009-4465-6005)). License: MIT (own code); third-party
files keep their own license (see `dados/THIRD_PARTY.md`).

> Code comments, logs and pre-registration files are in Portuguese; the technical note and this README are in
> English. Questions and issues are welcome.

## Context

[cyclosynth](https://github.com/mtweiden/cyclosynth) (Weiden et al., arXiv:2609.16659) synthesizes single-qubit
rotations over Clifford+√T with cost = n_T + 3·n_√T, about 24% cheaper than Clifford+T on the authors' targets.
The authors state that cost-optimality is not proved. This repository measures the gap.

## Results

### 1. Optimality gap of the published circuits (technical note: `NOTA_L3_L5.md`, PDF in `envio_autores/`)

Meet-in-the-middle exhaustive search over canonical Clifford+√T words (Forest et al., arXiv:1501.04944), exact
arithmetic in Z[ζ₁₆], on the authors' own published targets (`dados/u3.csv`):

| | ε = 1e-3 | ε = 1e-4 |
|---|---|---|
| Targets | 500 | 150 |
| A strictly cheaper circuit exists | 22 (4.4%) | 29 (19.3%) |
| Published circuit certified optimal | 478 | 120 (1 undetermined) |
| Excess over total published cost | 0.29% | 1.20% |

All 51 cheaper circuits are listed (`envio_autores/l3_*_testemunhas.csv`) and were re-verified with the authors'
unmodified cost function and 60-digit arithmetic (`l3_conferir.py`). Controls: the authors' distance column is
reproduced; the search never misses a published circuit; word counts match the counting recurrence.

### 2. Optimal mixed synthesis in Clifford+√T (`sonda_M_mistura.py`, `verifica_M50.py`)

Using the *complete* candidate list, a linear program picks the mixture of circuits with minimum expected cost
subject to a rigorous diamond-norm bound, (1/2)‖E − Id‖⋄ ≤ |Σ pᵢ cᵢ aᵢ| + Σ pᵢ |aᵢ|².

| Effective error | Feasible | Expected cost (mixture) | Published deterministic | Ratio |
|---|---|---|---|---|
| 1e-6 | 40/40 | 18.67 | 41.41 | **0.451** |
| 1e-7 | 20/20 | 23.70 | 50.60 | **0.468** |

- Verified in 50-digit arithmetic against an independent lower bound (ancilla-assisted state search):
  lower ≤ upper ≤ ε in 111/111 mixtures; where the linear term vanishes the bound is exact.
- Fair comparison against Clifford+T mixtures on the same targets (`sonda_M_soT.py`): cost ratio
  √T/T = 0.752 (1e-6) and 0.751 (1e-7). The √T advantage survives mixing, with no extra synergy.
- "Optimal" means optimal over the enumerated list (cost ≤ coverage, distance < search radius) under the
  linearized bound.

**Prior art.** Randomized synthesis is known: Campbell (2017), Hastings (2017), Akibue–Kato–Tani
(arXiv:2301.06307), and optimal Clifford+T mixtures by Morisaki–Sano–Akibue (arXiv:2510.05816). The contribution
here is the independent measurement in Clifford+√T and an open, verified benchmark.

### 3. A parameter-free slope estimate (`l5_lei_cambio.py`, exploratory)

From the counting model of Mooney, Hill and Hollenberg (Quantum 5, 396, 2021): cost ≈ (3/log₂x)·log₂(1/ε), with x
the root of Σ 2^(ℓ−2)·x^(−c_ℓ) = 1. It matches ten published slopes within 5% and gives 2.18·log₂(1/ε) for
Clifford+√T (fitted value in arXiv:2609.16659: 2.4).

## Reproducing

Windows or Linux, Python ≥ 3.12 with numpy, scipy, mpmath. The cyclosynth-dependent scripts (L4, mixing probes)
run in a Linux/WSL virtual environment with the `cyclosynth` package installed. Run from the repository root.

| Result | Command | Time (16 GB laptop) |
|---|---|---|
| Re-verify the 51 cheaper circuits | `python -X utf8 l3_conferir.py A B` | seconds |
| Exhaustive search | `python -X utf8 l3.py A` / `python -X utf8 l3.py B` | minutes to hours |
| L2 (denominator exponent vs. syllables) | `python -X utf8 l2.py` | minutes |
| Mixing probe | `python sonda_M_mistura.py A --n 40` / `B --n 20` | ~4 / ~3 min |
| Clifford+T mixing comparison | `python sonda_M_soT.py A --sqrt logs/sonda_M_A_20261120_042636.csv` | ~10 min |
| 50-digit verification | `python verifica_M50.py sonda_M_A_20261120_042636.csv sonda_M_B_20261120_043031.csv` | ~1 h |
| Slope estimate | `python -X utf8 l5_lei_cambio.py` | seconds |

Every run writes a timestamped log in `logs/`; the logs of the reported runs are included. Pre-registrations
(`PREREGISTRO_L*.md`) were written before the corresponding runs; later adjustments are recorded, with date and
reason, at the top of each script.

## Citation

If you use this code or data, please cite this repository and the cyclosynth paper (arXiv:2609.16659).
