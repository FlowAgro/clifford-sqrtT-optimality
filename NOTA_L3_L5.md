# Measured optimality gap of published Clifford+√T single-qubit circuits (technical note)

Author: Cristian Vitorino Rodrigues da Cunha (independent researcher, Uberlândia, MG, Brazil; ORCID 0009-0009-4465-6005) · 15 Nov 2026 · Code and data: this folder (`l2_ciclotomico/`)

## Summary

We ran an independent exhaustive search over Clifford+√T circuits for the Haar-random targets published
with cyclosynth (Weiden et al., arXiv:2609.16659; data file `u3.csv`). Cost is the authors' own:
cost = n_T + 3·n_√T.

- At ε = 1e-3, 22 of 500 targets (4.4%, 95% Wilson interval 2.9%–6.6%) admit a strictly cheaper circuit
  within the same tolerance; the other 478 published circuits are certified cost-optimal.
- At ε = 1e-4, 29 of 150 targets (19.3%, interval 13.8%–26.4%) admit a cheaper circuit; 120 are certified
  optimal and 1 is undetermined.
- The total excess cost is small: 0.29% at 1e-3 and 1.20% at 1e-4.
- All 51 cheaper circuits are listed explicitly and were re-verified through an independent code path.

The paper states that no cost-optimality is proved, so this note does not contradict it. It measures the
size of the gap and how it grows with precision.

## Method

Meet-in-the-middle exhaustive search. Every circuit is written as (canonical word L)·(canonical word R)·
(Clifford), with cost(L) ≤ a and cost(R) ≤ b; this covers all circuits of cost ≤ a + b − 2. Canonical words
follow Forest et al. (arXiv:1501.04944).

| Config | ε | Targets | a, b | Covers cost ≤ |
|---|---|---|---|---|
| A | 1e-3 | 0–499 | 12, 11 | 21 |
| B | 1e-4 | 0–149 | 16, 13 | 27 |

- Target: V = Rz(α)·Ry(β)·Rz(γ), as in the authors' `gather_u3.py`.
- Acceptance: d = sqrt(1 − |tr(U†V)|²/4) < ε, strict, no slack.
- Every circuit found is rebuilt as a gate sequence, its matrix is recomputed exactly in Z[ζ₁₆], and its
  cost and distance are recomputed from that matrix.
- The design was written down before the search code (`PREREGISTRO_L3.md`). The search has no tunable
  parameter that changes the outcome.

Controls, all passed in both configurations:

1. The authors' `distance` column is reproduced from their gates (relative error ≤ 1e-4; 0 failures).
2. For every target within reach, the search finds a circuit of cost ≤ the published cost (0 failures).
3. The number of generated words matches the counting recurrence (155,596 and 403,660).
4. The 24 Cliffords are distinct up to phase and have cost 0.

## Results

| | A (ε = 1e-3) | B (ε = 1e-4) |
|---|---|---|
| Targets | 500 | 150 |
| Cheaper circuit exists | 22 (4.4%) | 29 (19.3%) |
| Published circuit certified optimal | 478 | 120 |
| Undetermined | 0 | 1 (target 52) |
| Mean saving where suboptimal (T states) | 1.23 (max 3) | 1.62 (max 4) |
| Excess over total cost | 0.29% | 1.20% (47 of 3922) |

Independent verification (`l3_conferir.py`): cost computed with the authors' unmodified `_cost.py`;
distance recomputed with 60-digit arithmetic (mpmath), not sharing code with the search. 22/22 and 29/29
witnesses confirmed; largest d/ε is 0.998 (A) and 0.975 (B).

Witness files: `l3_A_testemunhas.csv`, `l3_B_testemunhas.csv` (ε, target index, published cost, new cost,
distance of the cheaper circuit, cheaper gates, published gates; no header row).

## Where the cheaper circuits are (exploratory, not pre-registered)

In all 51 witnesses the cheaper circuit has a larger denominator exponent than the published one (mean
+5.7 levels, range +1 to +14), fewer √T syllables (mean −2.7 in A, −3.0 in B) and more T syllables (+7.0,
+7.4). They are deeper and richer in T.

## Does a wider search window in cyclosynth recover them? (pre-registered, `PREREGISTRO_L4.md`)

cyclosynth at commit 0b1ffd0, called as the authors do (`sqrt_t=True, near_clifford_escape=False`).
Control: window 2 with no deadline reproduces the published cost in 149 of 150 targets at ε = 1e-4.

On the 29 suboptimal targets at ε = 1e-4 (published mean cost 27.79; optimal 26.17), no deadline:

| `lde_window` | At the optimum | Improved | Mean cost | Seconds per target |
|---|---|---|---|---|
| 4 | 8 (28%) | 9 | 27.17 | 7.9 |
| 6 | 9 (31%) | 11 | 27.07 | 16.3 |
| 8 | 12 (41%) | 16 | 26.86 | 26.0 |

Window 8 recovers 57% of the cost gap at about twelve times the run time of window 2 (26.0 s against
2.1 s per target). No returned circuit was
cheaper than our optimum, which cross-checks the exhaustive search.

At ε = 1e-7 (20 targets, optimum unknown, published mean 50.20), longer deadlines give 49.20 (window 2,
20 s), 48.90 (window 4, 20 s) and 48.70 (window 4, 60 s): up to 2.99% below the published cost, still
decreasing with time. With the default settings our own run gave 50.75, which is 1.10% above the published
cost (8 of 20 targets worse, 4 better), so the published runs are not reproduced exactly on our machine at
this precision. With n = 20 this is suggestive only.

Note: the source at this commit already contains a certified lower bound on cost from the denominator
exponent and an internal exhaustive certified mode not exposed in the Python interface. The gap measured
here is therefore a time-versus-cost trade-off in the default settings, not a missing capability.

## A counting estimate for the slope (exploratory)

Mooney, Hill and Hollenberg (Quantum 5, 396 (2021), arXiv:2005.05581, Appendix A, Eq. 20) count the
canonical sequences of cost ≤ C built from z-rotations of Clifford-hierarchy levels ℓ ≥ 3, with |T_ℓ| =
2^(ℓ−2) rotations of cost c_ℓ per level. They use this count to predict gate proportions. Its asymptotic
growth rate is x^C, where x is the root of

    Σ_ℓ 2^(ℓ−2) · x^(−c_ℓ) = 1.

Combined with the standard volume argument (about ε^(−3) distinct unitaries are needed), this gives the
typical-target estimate cost ≈ (3 / log2 x) · log2(1/ε), with no free parameter. We did not find this step
in that paper, where the slopes are reported as empirical fits; we regard it as a short corollary of their
model, not a new result.

Check against their fitted slopes (per log10(1/ε), their Tables 1 and 3; `l5_lei_cambio.py`):

| Gate set | Direct: predicted / measured | Magic states: predicted / measured |
|---|---|---|
| Set 1 | 9.97 / 10.46 ± 0.43 | 9.97 / 10.46 ± 0.43 |
| Set 2 | 6.62 / 6.89 ± 0.22 | 7.25 / 7.47 ± 0.15 |
| Set 3 | 5.86 / 6.05 ± 0.03 | 6.96 / 7.19 ± 0.12 |
| Set 4 | 5.40 / 5.33 ± 0.06 | 6.89 / 7.21 ± 0.25 |
| Set 5 | 5.02 / 4.84 ± 0.21 | 6.88 / 7.21 ± 0.25 |

All ten agree within 5%. In two cases (Sets 4 and 5, direct) the measured slope is below the estimate, by
1.3% and 3.5%, so the estimate should not be read as a strict bound on a fitted slope over a finite range.

For Clifford+√T with cost 3 per √T, the same equation is x³ = 2x² + 4, x = 2.594, giving 2.18·log2(1/ε),
against the 2.4·log2(1/ε) fitted in arXiv:2609.16659. Our optimal costs at 1e-3 and 1e-4 (18.71 and 25.83
on different target sets) correspond to about 2.14 per bit between the two precisions.

## Limits

- Exhaustive search does not reach ε ≤ 1e-5; the growth of the gap beyond 1e-4 is not tested.
- Configuration B covers 150 targets; one is undetermined.
- The counting section and the mechanism section are exploratory. The percentage reductions of Mooney et
  al. had been seen before the estimate was computed; the costs and slopes had not.
- Literature check for the slope estimate was limited to arXiv:2005.05581, 1501.04944, 2203.10064,
  2609.16659 and two web searches (15 Nov 2026). It may exist elsewhere.
- Run times were measured on a 12-thread laptop under WSL and are not comparable to the paper's.

## Files

`PREREGISTRO_L3.md`, `l3.py`, `l3_conferir.py`, `l3_A.log`, `l3_B.log`, `l3_A_testemunhas.csv`,
`l3_B_testemunhas.csv`, `PREREGISTRO_L4.md`, `l4_rodar.py`, `l4_analise.py`, `l4_A.csv`, `l4_B.csv`,
`l5_lei_cambio.py`, `l5_lei_cambio.log`.
