# Part VII — Geometric–thermodynamic certification

Model D6: compact box `Θ ⊂ ℝᵏ`, continuous `X: Θ → Conf_n`, continuous `U: Θ → ℝ`, `μ` with
continuous bounded density w.r.t. Lebesgue, `π(dθ) = Z⁻¹e^{−βU}μ(dθ)`, `A = {θ : P_{a,b}(X(θ))}`.
Nothing here concerns binding affinity; `Pr(A)` is an equilibrium population in a stated model.

## Lemma T0
`A = {θ : τ_{a,b}(X(θ)) ≥ 0}` is closed (F5 and S1: `τ` is 1-Lipschitz in `d_∞`), hence measurable.

## Theorem T1 (transfer principle: geometric certificates ⇒ probability bounds)
For **any** probability measure `P` on configurations and any region `U ⊆ Conf_n`:
if `U` is certified pocket-free then `P(pocket) ≤ P(X ∉ U)`; if every `X ∈ U` is certified to
have a pocket then `P(pocket) ≥ P(X ∈ U)`. For per-atom regions,
`P(X ∉ U(X⁰;ε,0)) ≤ Σᵢ P(‖xᵢ − x⁰ᵢ‖ > εᵢ)`. In the harmonic model (`ξᵢ ~ N(0,Σᵢᵢ)`),
`P(‖ξᵢ‖ > εᵢ) ≤ Q₃(εᵢ²/λ_max(Σᵢᵢ))` with the exact χ²₃ tail
`Q₃(s) = erfc(√(s/2)) + √(2s/π) e^{−s/2}`.
*Proof.* Inclusion of events; union bound; `‖ξᵢ‖² ≤ λ_max(Σᵢᵢ)‖zᵢ‖²` with `zᵢ ~ N(0,I₃)`. ∎
*Use.* Choosing non-uniform `εᵢ` (large for atoms whose motion is irrelevant, as certified by R4
and S3) makes the bound sharp in practice: a robust certificate converts directly into a
rigorous population bound under *any* stated probability model, with no sampling.

## Theorem T2 (branch-and-bound bounds on prevalence: soundness)
Let `{Q}` partition `Θ` into boxes with status `IN` (`Q ⊆ A`), `OUT` (`Q ∩ A = ∅`) or `UNK`, and
`U_lo(Q) ≤ U ≤ U_hi(Q)` on `Q`. Let `w_lo(Q) = μ(Q)e^{−βU_hi(Q)}`, `w_hi(Q) = μ(Q)e^{−βU_lo(Q)}`,
`Z_A^lo = Σ_IN w_lo`, `Z_A^hi = Σ_{IN∪UNK} w_hi`, `Z_N^lo = Σ_OUT w_lo`,
`Z_N^hi = Σ_{OUT∪UNK} w_hi`. Then
`L := Z_A^lo/(Z_A^lo + Z_N^hi) ≤ π(A) ≤ Z_A^hi/(Z_A^hi + Z_N^lo) =: Ū`.
*Proof.* `π(A) = Z_A/(Z_A+Z_N)` with `Z_A = ∫_A e^{−βU}dμ ∈ [Z_A^lo, Z_A^hi]`,
`Z_N ∈ [Z_N^lo, Z_N^hi]` (each `UNK` box is split between the two integrals, and both integrals
are bounded using the whole box). The map `(z_a, z_n) ↦ z_a/(z_a+z_n)` is increasing in `z_a` and
decreasing in `z_n`. ∎ (Lean: `PocketGenesis/Prevalence.lean`, `ratio_bounds`.)

The free-energy difference `ΔG = −β⁻¹ ln(Z_A/Z_N)` obeys
`−β⁻¹ln(Z_A^hi/Z_N^lo) ≤ ΔG ≤ −β⁻¹ ln(Z_A^lo/Z_N^hi)`.

## Theorem T3 (convergence under refinement)
Assume (i) `U` continuous and the energy enclosures have width → 0 with box diameter;
(ii) `μ(Z) = 0` for `Z = {θ : τ(X(θ)) = 0}`; (iii) the certification oracle is *asymptotically
complete*: for each `θ ∉ Z` there is `ρ(θ) > 0` such that every box of diameter `< ρ(θ)` containing
`θ` is certified `IN` (if `τ(θ) > 0`) or `OUT` (if `τ(θ) < 0`). If all `UNK` boxes are bisected
in every round, then `Ū_k − L_k → 0`.
*Proof.* The `UNK` mass at round `k` is `∫ 1[θ in an UNK box] dμ`; for `θ ∉ Z` the integrand is
eventually 0, so by dominated convergence and (ii) it tends to 0. By (i) and uniform continuity,
`Σ_{IN} (w_hi − w_lo) → 0` and similarly for `OUT`. Hence `Z_A^hi − Z_A^lo → 0`,
`Z_N^hi − Z_N^lo → 0`, and `L, Ū → π(A)` by continuity of the ratio (`Z > 0`). ∎
*Oracle completeness.* (iii) holds for the oracle built from R2/R5 (`IN`, in general position) and
Algorithm CL tests T1–T2 (`OUT`): if `τ(θ) < −3η`, every chamber candidate has `c < a − η` or an
escape path of clearance `> b + η`, which B2 + R1 certify robustly for configuration boxes of
`d_∞`-radius `< η`.
*Failure of (ii).* If the pocket exists only on a lower-dimensional set (e.g. exactly at a
tangency), `π(A) = 0` while `A ≠ ∅`; then `L = 0` always and `Ū → 0` — still correct.

## Proposition T4 (rate and cost)
If `θ ↦ X(θ)` is `ℓ_X`-Lipschitz into `d_∞`, `U` is `ℓ_U`-Lipschitz, boxes of `d_∞`-radius
`≤ |τ|/c₀` are always certified, and `μ{|τ∘X| ≤ s} ≤ κs`, then after refining to diameter `h`,
`Ū − L ≤ C (β ℓ_U h + κ c₀ ℓ_X h)` and adaptive refinement uses `O(κ c₀ ℓ_X h^{1−k})` boxes at the
finest level. Certified prevalence is therefore practical only for low-dimensional collective
coordinates (`k ≲ 3–4`); in high dimension the cost is exponential, matching the classical
lower bound `Ω(η^{−k})` for certified integration of Lipschitz functions (information-based
complexity). This is an honest limitation, not an artefact of the method.

## Theorem T5 (low-temperature limit: pocket free energy → minimal elastic work)
Harmonic model `ξ ~ N(0, k_BT·C)` on `range(H)`, `C = H⁺`, `X⁰` pocket-free. Let
`W* = inf{½ξᵀHξ : P_{a,b}(X⁰+ξ)}`. If `W*` equals the infimum over the interior of the pocket
set (regularity), then `lim_{T→0} k_BT ln π_T(A) = −W*`; hence
`ΔG_pocket(T) = −k_BT ln(π_T(A)/(1 − π_T(A))) → W*`.
*Proof.* Finite-dimensional Gaussian large deviations (Dembo–Zeitouni, Thm 2.2.30/Cramér; or
direct Laplace asymptotics) with rate `I(ξ) = ½ξᵀHξ`; `A` closed gives the upper bound, the
interior condition the lower bound. ∎ Combined with C6, certified `W*_lo ≤ W* ≤ W*_hi` yield
certified asymptotic bounds on the pocket free energy. The regularity hypothesis is necessary
(CE-T1).
