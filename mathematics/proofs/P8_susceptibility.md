# Part VIII — The pocket formation operator: derivatives, susceptibility, nonsmoothness

## Proposition K1 (critical values and their derivatives)
Let `s` be a critical point of `c_X` (i.e. `0 ∈ conv{uᵢ : i ∈ I}`, `uᵢ = (s−xᵢ)/‖s−xᵢ‖`,
`I = {i : fᵢ(s) = c(s)}`), equivalently `s ∈ relint conv{xᵢ : i ∈ I}` with
`Σ_{i∈I} λᵢuᵢ = 0`, `λᵢ > 0`, `Σλᵢ = 1`. Call `s` *nondegenerate* if the `uᵢ` are affinely
independent and the Hessian of `Σλᵢfᵢ` restricted to the tangent space of `{fᵢ = f_j ∀i,j∈I}` is
nonsingular. Then near `X` the critical point persists, its value `v(X) = c_X(s(X))` is
real-analytic, and
`∂v/∂xᵢ = −λᵢuᵢ`, `∂v/∂rᵢ = −λᵢ` for `i ∈ I`; all other partial derivatives vanish.
*Proof.* Implicit function theorem for the square system
`{fᵢ(s)−f_j(s) = 0 (i,j∈I), Σλᵢ∇fᵢ(s) = 0, Σλᵢ = 1}` in `(s,λ)`; nondegeneracy is exactly
invertibility of its Jacobian. Differentiating `v = fᵢ(s(X);X)` for each `i ∈ I` and averaging
with weights `λᵢ`: `dv = (Σλᵢ∇_yfᵢ)·ds + Σλᵢ ∂_Xfᵢ = Σλᵢ∂_Xfᵢ`, with `∂fᵢ/∂xᵢ = −uᵢ`,
`∂fᵢ/∂rᵢ = −1`. ∎
*Index.* `|I| = 1` cannot be critical; `|I| = 2` (midpoint-type) is index 1 for the distance
function, `|I| = 3` (gate, mouth window) index 2, `|I| = 4` (chamber centre) index 3
(Gershkovich–Rubinstein; Edelsbrunner's flow; Giesen–John flow complex). Chambers are index-3
points of `−(−c)`, i.e. local maxima of `c`; gates are index-2 saddles.

## Theorem K2 (pocket susceptibility on regular strata)
On the open set of configurations where a basin `W` has a nondegenerate chamber maximum `m`
(active set `I_m`, weights `μᵢ`) and a unique nondegenerate gate `s` (active set `I_s`, weights
`λᵢ`), the basin depth `D_W = c(m) − c(s)` is real-analytic with gradient
`∂D_W/∂xᵢ = −μᵢuᵢ^{(m)}·1[i∈I_m] + λᵢuᵢ^{(s)}·1[i∈I_s]`, `∂D_W/∂rᵢ = −μᵢ1[i∈I_m] + λᵢ1[i∈I_s]`.
The *susceptibility vector* `g = ∇_X D_W` has at most 7 nonzero atom blocks and
`Σᵢ ‖gᵢ‖ ≤ 2`, the infinitesimal form of S1's constant 2 (equality iff `I_m ∩ I_s = ∅`).
*Interpretation.* To deepen a pocket, push chamber atoms outward along `uᵢ^{(m)}` with weights
`μᵢ` and gate atoms inward toward the gate along `uᵢ^{(s)}` with weights `λᵢ`; every other atom
has zero first-order effect.

## Corollary K3 (linear response in the harmonic model)
With `Σ = k_BT·C`, to first order `Var(D_W) ≈ gᵀΣg = k_BT Σ_k (g·v_k)²/ω_k²` (normal modes
`v_k`, frequencies `ω_k`). The least-work deformation deepening the pocket by `Δ` within the
stratum is `ξ* = Δ·Cg/(gᵀCg)` with work `Δ²/(2gᵀCg)`. This is a *heuristic estimate* of `W*`;
it becomes a *rigorous upper bound* only after `X⁰ + ξ*` is certified (C6).

## Proposition K4 (nonsmoothness and its structure)
`D_X(y)` and `τ_{a,b}` are globally Lipschitz (S1), so Clarke generalized gradients exist
everywhere. At a *gate switch* (two gates of `W` at equal height) `e_W` is the maximum of two
analytic functions and `D_W` the minimum of two, with Clarke gradient the convex hull of the
two-sided gradients; at *births/deaths* (a chamber/gate pair coalescing) the basin itself
disappears and the operator must be read in the persistence-diagram metric (Part IV), in which
it remains 1-Lipschitz.
*Limitation (counterexample CE-K1).* The softest normal mode need not open any pocket: a mode
that translates the pocket region rigidly has `g·v = 0` (translation invariance), so spectral
softness alone does not predict pocket formation.
