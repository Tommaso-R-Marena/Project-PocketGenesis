# Part II — Stability: sharp global bounds, local bounds, and what is *not* stable

Throughout, `X, Y ∈ Conf_n` share atom labels, `εᵢ = ‖xᵢ−yᵢ‖ + |rᵢ−r′ᵢ|`, `ε = maxᵢ εᵢ`, and both
use a common bulk `B_R` with `R ≥ max(ρ_X, ρ_Y)` (legitimate by F1(b)).

---

## Theorem S1 (global stability with optimal constants)
`‖c_X − c_Y‖_∞ ≤ ε`, `‖E_X − E_Y‖_∞ ≤ ε`, `‖D_X − D_Y‖_∞ ≤ 2ε`, `|τ_{a,b}(X) − τ_{a,b}(Y)| ≤ ε`.
The constants `1, 1, 2` are attained (Example S1 attains all three simultaneously).

*Proof.* `|f_i^X(z) − f_i^Y(z)| ≤ ‖xᵢ−yᵢ‖ + |rᵢ−r′ᵢ| ≤ ε` for every `z`, and `min` is 1-Lipschitz
for the sup-norm, so `|c_X(z) − c_Y(z)| ≤ ε`. For a fixed `γ ∈ Γ(y)` (the same family for both
configurations because the bulk is common), `|min_s c_X(γ(s)) − min_s c_Y(γ(s))| ≤ ε`; a supremum
of a family of functions that are pointwise within `ε` is within `ε`. `D` by the triangle
inequality; `τ` because `y ↦ min(c−a, b−E)` moves by at most `ε` pointwise. ∎

### Example S1 (sharpness — "octahedral cage with a chamber-control atom")
`X`: six cage atoms at `±ρe_k` (`k=1,2,3`) of radius `r`, plus one atom `q` at `−d u`,
`u = (1,1,1)/√3`, of radius `r_s`; numerically `ρ=10, r=1.5, d=9, r_s=1.5`; point `y=0`.
`Y`: cage radii `r+ε`, `q` moved to `−(d+ε)u`; `0 < ε ≤ ½`. Then `d_∞(X,Y) = ε` and

| quantity | `X` | `Y` | change |
|---|---|---|---|
| `c(0)` | `d − r_s = 7.5` | `7.5 + ε` | `+ε` |
| `E(0)` | `g := ρ√(2/3) − r ≈ 6.66497` | `g − ε` | `−ε` |
| `D(0)` | `7.5 − g` | `7.5 − g + 2ε` | `+2ε` |

*Verification.* `c_X(0) = min(ρ−r, d−r_s) = 7.5`; `c_Y(0) = min(8.5−ε, 7.5+ε) = 7.5+ε`.
*Upper bound on E:* `K = ∂conv{±ρe_k}` separates `0` from the bulk. On the face with vertices
`ρe₁,ρe₂,ρe₃` the three positive vertices are the nearest cage centres
(`‖z+ρe₁‖² − ‖z−ρe₁‖² = 4ρz₁ ≥ 0` there), and the maximum over an equilateral triangle of the
distance to its nearest vertex is the circumradius `ρ√2/√3`, attained at the centroid. Hence
`max_K c_cage,X = g`, and `max_K c_Y ≤ max_K c_cage,Y = g − ε` (every cage term drops by exactly
`ε`; extra atoms only lower `c`). By duality (Theorem B1) `E_X(0) ≤ g`, `E_Y(0) ≤ g−ε`.
*Lower bound on E:* along the ray `tu`, `‖tu−ρe₁‖² = t² − 2ρt/√3 + ρ²` has minimum `2ρ²/3`, the
negative vertices are farther, and the `q`-term is `t + d − r_s ≥ 7.5 ≥ g`; so `E_X(0) ≥ g`.
For `Y`, S1 gives `E_Y(0) ≥ E_X(0) − ε = g − ε`. ∎ (This example is a regression test:
`tests/test_sharpness.py`.)

## Theorem S2 (nonuniform local stability of the clearance)
Let `m(z) = minᵢ (f_i^X(z) + εᵢ)` and `N(z) = { i : f_i^X(z) − εᵢ ≤ m(z) }` (atoms that can be
nearest to `z` under the perturbation). Then
`c_X(z) − max_{i ∈ N(z)} εᵢ ≤ c_Y(z) ≤ c_X(z) + min_{j ∈ I_X(z)} ε_j`,
where `I_X(z) = argminᵢ f_i^X(z)`. Write `ε̄(z) := max_{i∈N(z)} εᵢ`.

*Proof.* Upper: `c_Y(z) ≤ f_j^Y(z) ≤ f_j^X(z) + ε_j = c_X(z) + ε_j` for `j ∈ I_X(z)`. Lower: for
`i ∉ N(z)`, `f_i^Y(z) ≥ f_i^X(z) − εᵢ > m(z) ≥ c_Y(z)`, so the minimum defining `c_Y(z)` is attained
in `N(z)`, where `f_i^Y ≥ f_i^X − εᵢ ≥ c_X(z) − ε̄(z)`. ∎

## Theorem S3 (local stability of escape clearance — only path and barrier atoms matter)
Let `γ*` be an optimal escape path of `y` for `X` (F3) and `K*` a compact separating set with
`max_{K*} c_X = E_X(y)` (Theorem B1). Then
`E_X(y) − max_{z∈γ*} ε̄(z) ≤ E_Y(y) ≤ E_X(y) + max_{z∈K*} min_{j∈I_X(z)} ε_j`.
In particular `E_Y(y) = E_X(y)` whenever all atoms that can be nearest to a point of `γ*` or of
`K*` are fixed; and the pocket at a chamber witness `y` persists
(`P_{a,b}(Y)` with the same witness) whenever
`max_{i ∈ N(y)} εᵢ ≤ c_X(y) − a` and `max_{z∈K*} min_{j∈I_X(z)} ε_j ≤ b − E_X(y)`.

*Proof.* `E_Y(y) ≥ min_s c_Y(γ*(s)) ≥ min_s (c_X(γ*(s)) − ε̄(γ*(s)))`. `K*` separates `y` from the
(common) bulk in both configurations (separation is a property of the set, not of `X`), so by B1
`E_Y(y) ≤ max_{K*} c_Y` and S2's upper bound applies pointwise. ∎

*Why this matters.* Global bounds charge every atom's displacement to every pocket. S3 charges only
the "pocket-control atoms" `N(γ*) ∪ N(K*) ∪ N(y)`. Example: a hinge motion that moves 90% of a
protein by 5 Å while leaving the atoms near a pocket's chamber, escape path and barrier fixed
changes nothing (global bound: 5 Å; local bound: 0). See `COUNTEREXAMPLES.md`, CE-S2.

## Theorem S4 (two-parameter interleaving and certified volume bounds)
For `s > t` define the **chamber–mouth bifiltration** `P_X(s,t) = { y : c_X(y) > s, E_X(y) ≤ t }`
(bounded: by F6 it lies in the union of basins ⊆ `B(o,R)`). It is decreasing in `s`, increasing in
`t`, and
`P_X(s+ε, t−ε) ⊆ P_Y(s,t) ⊆ P_X(s−ε, t+ε)` whenever `s − t > 2ε`.
Hence `v_X(s,t) = Leb(P_X(s,t))` satisfies `v_X(s+ε,t−ε) ≤ v_Y(s,t) ≤ v_X(s−ε,t+ε)`.

*Proof.* If `c_X(y) > s+ε` then `c_Y(y) > s`; if `E_X(y) ≤ t−ε` then `E_Y(y) ≤ t`. The other
inclusion is symmetric. `c, E` are continuous so the sets are Borel; monotonicity of Lebesgue
measure gives the volume bounds. ∎

## Counterexample CE-S1 (single-scale pocket volume is not stable — not even continuous)
The cavity volume at one probe radius, `b ↦ v(b,b) = Leb(bounded components of {c > b})`, is a
discontinuous function of the configuration. In the octahedral cage (no extra atom; `ρ=10,
r=1.5`) with `b = g = ρ√(2/3) − r`, `v_X(b,b) ≥ (4π/3)(ρ − r − g)³ ≈ 25.9` (a ball around the
centre is enclosed), while shrinking all radii by any `ε > 0` gives `g_Y = g + ε > b`, the interior
connects to bulk through the eight windows, and `v_Y(b,b) = 0`. Theorem S4 is the correct
replacement: the stable object is the two-parameter function `v(s,t)` under the shift
`(s,t) ↦ (s±ε, t∓ε)` — exactly analogous to why persistence diagrams, not Betti numbers, are
stable.

## Proposition S5 (the pocket formation operator is 1-Lipschitz into diagram space)
The map `X ↦ Dgm(X)` (relative pocket diagram, Part IV) satisfies
`d_B(Dgm(X), Dgm(Y)) ≤ d_∞(X,Y)`, is equivariant under isometries and invariant under
relabelling. The constant 1 is optimal (single-chamber cage scaled radially).
*Proof.* Part IV, Theorem H1(c), with S1. Optimality: a cage whose chamber value and gate value
both depend on all atoms shifts the single off-diagonal point by exactly `ε` in both coordinates
under a uniform radius change. ∎
