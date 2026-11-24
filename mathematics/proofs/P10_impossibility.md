# Part X — Fundamental limitations

## Theorem I1 (accessibility is not local; certificates must be as large as the pocket)
(a) For every `R > 0` there are configurations `X′, X″` that coincide on `B(y,R)` such that `y`
is a chamber witness for `X′` but not for `X″`. (b) If `K ⊆ {c ≤ b}` separates `y` from the bulk
and `W` is the component of `{c > b}` containing `y`, then `W` lies in the bounded region enclosed
by `K`, so `diam K ≥ diam W`; a barrier certificate for a pocket whose accessible region has
diameter `L` involves atoms spread over diameter `≥ L − 2(b + r_max)`.
*Proof.* (a) A straight tunnel of length `> R` from the chamber to the surface, closed at the far
end by one extra ball in `X′` and open in `X″`. (b) `W` is connected and disjoint from `K`. ∎
*Reconciliation with R4.* Confinement is local *to the barrier atoms*, but the barrier is as large
as the confined region; locality cannot be measured from the chamber.

## Theorem I2 (binary pocket classifiers cannot be stable; abstention is forced)
`Conf_n` is connected, so every continuous map `Conf_n → {0,1}` is constant. If a procedure
receives `X̃` with the guarantee `d_∞(X̃, X) ≤ ε` and must answer soundly, it must abstain on the
*ambiguous set* `Amb_ε = {X̃ : ∃X⁺,X⁻ ∈ B_ε(X̃), P(X⁺), ¬P(X⁻)}`, and `Amb_ε ⊆ {|τ(X̃)| ≤ ε}` (S1).
Our checker abstains (`UNKNOWN`) only inside `{|τ| ≤ ε + δ_cert}` where `δ_cert` is the
certificate resolution (R5, B3), i.e. it is optimal up to `δ_cert`. ∎

## Theorem I3 (finite samples cannot certify absence or rarity)
(a) Without a measure or dynamical model, no finite set of observed conformations implies
`¬P_{a,b}` on `M` (any unobserved `X ∈ M` may carry the pocket). (b) With `N` i.i.d. samples from
`π` and no pocket observed, the exact one-sided `1−α` Clopper–Pearson bound is
`π(A) ≤ 1 − α^{1/N}` (`≈ 3/N` for `α = 0.05`); `π(A) = 0` is never implied, and the statement is a
confidence statement, not a certificate. Certificates of rarity require Part VII (a model and an
exhaustive certified decomposition) or T1 (a model and a certified region).
*Proof.* (a) trivial; (b) the binomial likelihood of zero successes. ∎

## Theorem I4 (geometry does not determine prevalence)
For any measurable `A ⊂ Θ` with `0 < μ(A) < μ(Θ)` and any `p ∈ (0,1)` there is a continuous
energy `U` with `π_U(A) = p`. Hence equilibrium pocket populations cannot be inferred from the
configuration space and the pocket predicate alone. ∎ (Take `U = λ·dist(·, A_δ)` type functions
and use continuity of `λ ↦ π_{λ}(A)` with limits `μ`-proportions at `λ = 0` and `1` as `λ → ∞`
after smoothing; intermediate value theorem.)

## Theorem I5 (geometric accessibility is neither necessary nor sufficient for energetic accessibility)
With entry barriers `B_V` (D7): for any configuration and any `y`, there are continuous probe
potentials `V` with `B_V(y)` arbitrarily large although `E_X(y) > ρ` (soft long-range repulsion
along all paths), and potentials with `B_V(y)` arbitrarily small although `E_X(y) < ρ`
(hard-sphere radii exceeding effective softness). Only within a stated potential is accessibility
defined. ∎

## Theorem I6 (static endpoints do not determine dynamics)
There are trajectories `X₁, X₂` with `X₁(0) = X₂(0)`, `X₁(T) = X₂(T)`, both pocket-free at the
endpoints, such that a pocket forms and closes along `X₁` and never along `X₂`. ∎ (Close a gate
transiently.) Endpoint comparison (apo/holo) can therefore never certify the absence of
intermediate events; Dy2 needs samples with Lipschitz control.

## Pointers to impossibility results proved elsewhere
- No global continuous pocket identity: Theorem M1 (Part V).
- No exact single-radius criterion for rigid probes: Theorem Q2 (Part IX).
- Discrete-state cryptic pocket search is NP-complete: Theorem C5 (Part VI).
- Infinitely many events for smooth (non-tame) motions: Theorem Dy4(b).
- Single-scale pocket volume is discontinuous: CE-S1 (Part II).
- The stability lower bound on formation displacement can be infinitely weak: Theorem C3.
