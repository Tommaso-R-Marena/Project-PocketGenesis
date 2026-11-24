# Part VI — Cryptic pockets: obstructions, exact enclosure numbers, hardness, elastic work

`δ*(X⁰; M) = inf{ d(X⁰,X) : X ∈ M, P_{a,b}(X) }` (D3.9). Unless stated otherwise `M` is the
*independent-atom model*: `M = { X : ‖xᵢ − x⁰ᵢ‖ ≤ δᵢ }` with radii fixed, or `M = Conf_n` with the
metric `d_∞`. Every lower bound proved for a family `M` holds for every subfamily (sterically
feasible, energy-bounded, ...): `M′ ⊆ M ⇒ δ*(M′) ≥ δ*(M)`.

## Theorem C1 (stability lower bound — calibration)
`δ*_∞(X⁰) ≥ −τ_{a,b}(X⁰)`; site-specifically, a pocket with witness `y` in some `X` forces
`d_∞(X⁰,X) ≥ δ_stab(y) := max(a − c_{X⁰}(y), E_{X⁰}(y) − b)`.
*Proof.* S1. ∎

The stability bound treats the two requirements of a pocket — *clear the chamber* and *close the
mouth* — independently. The next results exploit that the same atoms must do both.

## Theorem C2 (chamber-constrained escape-ray obstruction)
Let `Y` be a box with centre `y_c` and half-diagonal `ρ_Y`, per-atom budgets `δᵢ`,
`d_i = ‖x⁰ᵢ − y_c‖`, `uᵢ = (x⁰ᵢ − y_c)/dᵢ`, and
`devᵢ = arcsin((δᵢ+ρ_Y)/dᵢ)` if `dᵢ > δᵢ + ρ_Y` (else `π`),
`κᵢ = arcsin( (rᵢ + b) / max(a + rᵢ, dᵢ − δᵢ − ρ_Y) )`, `Θᵢ = min(π, devᵢ + κᵢ)`.
If some unit vector `u*` satisfies `∠(u*, uᵢ) > Θᵢ` for every `i`, then **no** `X` with
`‖xᵢ − x⁰ᵢ‖ ≤ δᵢ` has a chamber witness in `Y`.
*Proof.* Suppose `X` has a witness `y ∈ Y`. The relative position `xᵢ − y` differs from
`x⁰ᵢ − y_c` by at most `δᵢ + ρ_Y`, so (tangent-cone bound) `∠(xᵢ − y, uᵢ) ≤ devᵢ`. The chamber
condition gives `‖xᵢ − y‖ ≥ a + rᵢ`, and the triangle inequality `‖xᵢ − y‖ ≥ dᵢ − δᵢ − ρ_Y`; the
"shadow" of `B̄(xᵢ, rᵢ+b)` seen from `y` is the cap of angular radius
`arcsin((rᵢ+b)/‖xᵢ−y‖) ≤ κᵢ` about `(xᵢ−y)/‖xᵢ−y‖`. Hence the shadow lies in `Cap(uᵢ, Θᵢ)`, the
ray `y + ℝ_{≥0}u*` misses every inflated ball, `c_X > b` along it (and `c_X(y) ≥ a > b`), the ray
reaches the bulk, and `E_X(y) > b` — contradiction. ∎
*Remark.* The test is verified with rigorous ball arithmetic (Arb via `python-flint`), not
floating point. Ordinary robust escape paths (R1b with inflation `δᵢ + ρ_Y`) are a second sound
test; C2 is strictly stronger exactly when atoms that must vacate the chamber are the ones that
would have to close the mouth (Theorem C3).

## Theorem C4 (enclosure number = spherical covering number)
For `a > b`, `r ≥ 0`, `r + b > 0`, let `θ* = arcsin((r+b)/(a+r))` and let `N_cov(θ)` be the least
number of closed spherical caps of angular radius `θ` covering `S²`. Then
**the least number of balls of radius `r` forming a configuration satisfying `P_{a,b}` is exactly
`N_cov(θ*)`.** If radii vary in `[0, r_max]`, the least number is `N_cov(θ*(r_max))`.

*Proof.* (≥) Let `y` be a chamber witness. Every ray from `y` lies in `Γ(y)` after truncation, so
by F3 it meets `U_b`. Each ball has `‖xᵢ − y‖ ≥ a + rᵢ > rᵢ + b`, so its shadow is a cap of radius
`arcsin((rᵢ+b)/‖xᵢ−y‖) ≤ arcsin((rᵢ+b)/(a+rᵢ)) ≤ θ*(r_max)` (the map `r ↦ (r+b)/(a+r)` is increasing
because `a > b`). The shadows cover `S²`.
(≤) Let `p₁,…,p_N` be centres of an optimal cover by caps of radius `θ*`. Put balls of radius `r`
at `d·p_k`, `d = a + r`, and let `R = r + b`, `ρ′ = √(d² − R²) = d cos θ*`. For a unit vector `v`,
`‖ρ′v − d p_k‖² ≤ R² ⇔ v·p_k ≥ (ρ′² + d² − R²)/(2ρ′d) = ρ′/d = cos θ*`. So `B̄(d p_k, R)` meets the
sphere `S(0,ρ′)` exactly in the cap of radius `θ*` about `p_k`, the balls cover `S(0,ρ′)`, which
separates `0` from the bulk, hence `E(0) ≤ b` by B1; and `c(0) = d − r = a`. ∎

*Explicit bounds.* Fejes Tóth (1943): `N` caps of radius `θ` can cover `S²` only if
`cos θ ≤ cot(Nπ / (6(N−2))) / √3`, sharp for `N = 4, 6, 12` (tetrahedron 70.53°, octahedron
54.74°, icosahedron 37.38°); also `N ≥ 2/(1 − cos θ)` (area). Representative values for carbon-like
radius `r = 1.7` Å (computed in `experiments/synthetic/enclosure_numbers.py`):

| chamber `a` (Å) | mouth `b` (Å) | `θ*` | area bound | Fejes Tóth bound on # atoms |
|---|---|---|---|---|
| 1.6 | 1.4 | 69.95° | 4 | 5 |
| 3.0 | 1.4 | 41.27° | 9 | 10 |
| 4.0 | 1.4 | 32.95° | 13 | 16 |
| 4.0 | 0.5 | 22.70° | 26 | 32 |
| 6.0 | 0.5 | 16.60° | 48 | 59 |

(First row: `69.95°` is just below the tetrahedral covering radius `70.53°`, so four atoms of
radius 1.7 Å cannot confine a 1.6 Å chamber against a water-sized mouth — five are needed.)

## Theorem C3 (the stability bound can be arbitrarily weak)
There are configurations with `δ_stab ≤ 1` (Theorem C1 predicts a pocket might form after moving
atoms by 1 Å) but `δ* = ∞` (no motion of the atoms ever forms the pocket).
*Example.* Twelve balls of radius `1.5` at the vertices of a regular icosahedron of circumradius
`10` about `0`. Here `c(0) = 8.5` and the gate value is `g = 10 sin(37.377°) − 1.5 ≈ 4.5706`
(face-circumcentre clearance; Example S1's argument). For `a = 9.5`, `b = g − 1 ≈ 3.5706`,
`δ_stab(0) = 1`. But `θ* = arcsin(5.0706/11) ≈ 27.45° < 37.38°`, and Fejes Tóth gives
`N_cov(θ*) ≥ 22 > 12`: by C4 no configuration of twelve such balls satisfies `P_{a,b}`. ∎
Adding atoms at distance `λ + (a+r)` makes `δ*` finite but `≥ λ` while `δ_stab` stays `≤ 1`.

## Theorem C5 (discrete-state pocket formation is strongly NP-complete)
*Problem DSPF.* Input: rationals `a > b ≥ 0`; fixed balls `F`; groups `G₁,…,G_K`, each with two
states, each state a list of balls (a rotamer-library-like model). Question: is there a choice of
states such that `F ∪ (chosen balls)` satisfies `P_{a,b}`?
**DSPF is NP-complete; it is NP-hard even with all radii 1, integer coordinates bounded by
`O(√size)`, non-interacting groups (no steric constraints), two states per group, and under the
promise that either `τ_{a,b} ≥ 0.08` or `τ_{a,b} ≤ −0.5`.**

*Membership.* For fixed states, `P_{a,b}` is decided in polynomial time: compute the regular
triangulation of `(xᵢ,(rᵢ+b)²)` and the alpha complex at level 0 exactly over ℚ, extract the
bounded components `Q` of `ℝ³∖U_b` as unions of open cells (Theorem B5, step 3), and decide
`∃y ∈ Q : ‖y−xᵢ‖² ≥ (a+rᵢ)² ∀i` — an existential first-order sentence in a fixed number (3) of
variables, decidable in polynomial time (Renegar 1992; Basu–Pollack–Roy).
*Hardness (from 3-SAT).* Let `M` be even with `M ≥ 4√m + 16`. **Shell:** unit balls at all points
of `(2ℤ)³` on the boundary of the cube `C = [0,2M]³`, except `m` *window* points `w_c` on the face
`z = 0` with `x, y ∈ [8, 2M−8]` and pairwise distance `≥ 8`. **Groups:** for each variable `x_j`,
one ball per occurrence; in state `TRUE` the balls of positive occurrences sit at the windows of
their clauses and the others at (distinct) existing shell lattice points (where they change
nothing); in state `FALSE`, vice versa. Set `b = ½`, `a = M − 2`.
(⇐) If every window is plugged, the balls of radius `1.5` centred at the full lattice cover `∂C`
(covering radius of the square lattice of spacing 2 is `√2 < 1.5`), so `∂C ⊆ U_b` separates the
centre `y = (M,M,M)` from the bulk: `E(y) ≤ max_{∂C} c ≤ √2 − 1 < b − 0.08`, and
`c(y) = M − 1 = a + 1`.
(⇒) If window `w` is open, take any `z` with `c(z) ≥ a`. *Outside `C`:* all centres lie in the
convex set `C`, so moving along the outward normal of `C` at `z` increases every distance;
`E(z) = c(z) ≥ a > b`. *Inside `C`:* `c(z) ≥ M − 2` forces distance `≥ √((M−1)² − 2)` from every
face, so `z` is within `1.1` of the centre in each coordinate (for `M ≥ 16`); the path `z → w + (0,0,M) → w − (0,0,2) → −∞·e₃` keeps
clearance `≥ 1`: on the horizontal leg the nearest centres are at distance `≥ M − 2` (bottom/top)
and `≥ 7` (sides); on the vertical leg through the empty lattice point the nearest centres are the
four lattice neighbours at distance `√(4 + h²) ≥ 2`, other windows' plugs are `≥ 8` away in-plane;
below the cube distances increase. So `E(z) ≥ 1 = b + ½`, and `τ_{a,b} ≤ −½`.
Thus `P_{a,b}` holds iff every clause has a true literal. Sizes and coordinates are polynomial
(in fact `O(m)` balls with coordinates `O(√m)`). ∎

*Consequences and limits.* Exhaustive search over discrete side-chain libraries for cryptic
pockets cannot be polynomial unless P = NP, and the hardness survives any additive numerical
tolerance below `0.08` in the margin. **Limitation:** variable groups in the reduction are
spatially non-local (a group touches all windows of its variable's clauses). Hardness for *local*
groups with steric coupling (planar 3-SAT on the shell surface with wire gadgets) is Conjecture
CJ-3. The continuous independent-atom problem ("is `δ*_∞ ≤ δ`?") is open (Conjecture CJ-2).

## Theorem C6 (certified bounds on minimal elastic work)
Harmonic model `X = X⁰ + ξ`, `ξ ∈ range(H)`, energy `½ξᵀHξ`, `C = H⁺`. Let
`εᵢ(W) = √(2W λ_max(C_{ii}))` (`C_{ii}` the 3×3 diagonal block). If the region `U(X⁰; ε(W), 0)`
is certified pocket-free (all chamber candidates covered by T1–T3 of Algorithm CL), then
`W* ≥ W`. If `X⁰ + ξ̂` carries a valid PRESENT certificate, then `W* ≤ ½ξ̂ᵀHξ̂`.
*Proof.* `max{ u·ξᵢ : ½ξᵀHξ ≤ W, ξ ∈ range H } = √(2W uᵀC_{ii}u)` (Cauchy–Schwarz in the
`H`-inner product), so the energy ellipsoid lies in the per-atom region. ∎

## Algorithm CL (certified lower bound on formation displacement)
Input `X⁰, a, b`, budget profile `δᵢ = δ·wᵢ`. Search region `B̄(o, ρ_X + δ_max)` (beyond it
`E = c` by F1). Recursively subdivide boxes `Y`; a box is *cleared* by
(T1) chamber infeasible: `∃i : ‖x⁰ᵢ − y_c‖ + ρ_Y + δᵢ < a + rᵢ` (exact: compare squares);
(T2) robust escape: `B̄(y_c,ρ_Y)` clear at level `b` with inflation `δᵢ`, plus a robust path
certificate from `y_c` (R1b with inflation `δᵢ + ρ_Y`... only `δᵢ` needed after joining);
(T3) Theorem C2 with a candidate set of directions.
If all leaves are cleared, return "`δ*(profile w) > δ`" (PROVED); else UNKNOWN with the uncleared
boxes. Bisection on `δ` gives the best certified lower bound. Soundness: C1–C2, R1, B1.
