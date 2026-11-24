# Part IX — Beyond spherical probes

A rigid probe is a compact `K ⊂ ℝ³`; a placement `(p,R) ∈ SE(3)` occupies `p + RK`. It is *free*
if `p + RK` avoids the closed atoms; *accessible* if joined to a free placement in the bulk by a
path of free placements. Let `B̄(z₀,ρ_in) ⊆ K ⊆ B̄(z₁,ρ_out)`.

## Theorem Q1 (sandwich: spherical probes give rigorous one-sided tests)
(a) If `E_X(p + Rz₁) > ρ_out` then `K` at `(p,R)` is accessible (translate along an escape path of
the circumscribed centre at fixed orientation).
(b) If `K` at `(p,R)` is accessible then `E_X(p + Rz₀) > ρ_in` (follow the inscribed ball's centre
along any free path; compactness gives a positive margin).
*Proof.* (a) `p(s) + RK ⊆ B̄(p(s)+Rz₁, ρ_out)` and `c > ρ_out` there means the ball avoids the
closed atoms. (b) along a free path `(p(s),R(s))`, `B̄(q(s),ρ_in) ⊆` probe with
`q(s) = p(s)+R(s)z₀`, so `c(q(s)) > ρ_in`; the minimum over the compact path is `> ρ_in`, and
`q(1)` is in the bulk after enlarging `R` (F1). ∎

## Theorem Q2 (no single-radius criterion is exact for non-spherical probes)
Let `K` be a solid rod (cylinder of radius `w`, length `L > 6w`) with reference point `z` its
centre. For every `ρ ∈ ℝ` there is a union of balls for which "`E_X(p+Rz) > ρ`" disagrees with
accessibility of `K` at `(p,R)`. (Other reference points need only minor changes of the
constructions; we state the natural case.)
*Proof.* *Case `ρ > w`.* Enclose the rod in a closed shell whose only exit is a straight channel,
aligned with the rod axis, with bottleneck clearance `w′ ∈ (w, ρ]`. The rod slides out axially
(accessible), but every escape path of any point crosses the narrowest cross-section, where
`c ≤ w′`; so `E(p+Rz) ≤ w′ ≤ ρ` and the criterion says "inaccessible".
*Case `ρ ≤ w`.* Use instead an L-shaped corridor (90° bend) with medial clearance `w + η`. A ball of
radius `w + η/2` follows the medial curve around the bend, so `E(p+Rz) ≥ w + η/2 > ρ` for `z`
near the axis (criterion: "accessible"). But the rod's axis must stay in the thin tube of radius
`η` about the medial curve, and a rigid segment of length `L` cannot turn a 90° bend in a tube of
width `2η` once `L > 2√2·2η` (ladder-around-a-corner bound); so for `η < L/6` the rod is trapped. ∎
*Consequence.* Orientation-dependent passage is a genuine SE(3) phenomenon; spherical
calculations give only the bracket `[ρ_in, ρ_out]`.

## Proposition Q3 (translating convex probes: the whole theory transfers)
For convex `K` with `0 ∈ int K` at fixed orientation, define the gauge clearance
`c_K(y) = sup{t ≥ 0 : (y + tK) ∩ Ω(X) = ∅}` (with `c_K ≤ 0` inside). Then `c_K` is continuous,
the translation-accessibility functional `E_K` (max–min of `c_K` over escape paths) satisfies the
duality B1, the basin decomposition F6, path and barrier certificates (cells contained in
`Ω ⊕ (−tK)` replace cells contained in inflated balls), and stability
`‖c_{K,X} − c_{K,Y}‖_∞ ≤ ε/ρ_in(K)` in units of `K`-scaling.
*Proof.* Every proof in Parts I–III uses only continuity of `c`, path concatenation and
compactness (B1–B3), or Lipschitz bounds; the Euclidean ball of radius `ε` lies in `(ε/ρ_in)K`. ∎
(The Čech-cycle certificate does not transfer verbatim: B4a is specific to Euclidean balls.)
