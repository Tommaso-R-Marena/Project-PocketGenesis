# Part III — Escape duality, exact certificates, robust and intrinsic certificates

## Theorem B1 (escape duality: max–min paths = min–max separating sets)
Call a compact `K ⊂ ℝ³` *separating for y* if every `γ ∈ Γ(y)` meets `K`. Then
`E_X(y) = min { max_{z ∈ K} c_X(z) : K separating for y }`, and the minimum is attained.
The same holds for any continuous potential `V` (D7): `B_V(y) = max_K min_{z∈K} V(z)`.

*Proof.* (≤) If `K` separates and `γ ∈ Γ(y)`, then `γ` meets `K` at some `z`, so
`min_γ c ≤ c(z) ≤ max_K c`; take `sup_γ`. (≥, attainment) Let `e = E(y)`. If `c(y) = e`, `K = {y}`
works. Otherwise `c(y) > e`; let `W` be the component of the open set `{c > e}` containing `y`.
As in F6(2), `W ∩ B_R = ∅` (else `E(y) > e`), so `W` is bounded and `K := ∂W` is compact. Every
`γ ∈ Γ(y)` leaves `W`: the closed set `{s : γ(s) ∉ W}` is nonempty (it contains `1`), its minimum
`s₁` is `> 0`, and `γ(s₁) ∈ W̄ ∖ W = ∂W`. Points of `∂W` are not in `{c > e}` (a boundary point of
a component of an open set lies outside the open set), so `max_{∂W} c ≤ e`. ∎

The discrete certificates below are finite, exactly checkable realizations of the two sides of B1.

## Theorem B2 (polygonal path certificates: soundness and completeness)
*Soundness.* If `p₀ = y, …, p_m` satisfy `‖p_m−o‖ ≥ R` and `dist(xᵢ,[p_k,p_{k+1}]) ≥ rᵢ + b′` for
all `i,k` (where `rᵢ + b′ ≥ 0`), then `E_X(y) ≥ b′`.
*Completeness.* If `E_X(y) > b′` there is such a certificate with rational vertices.

*Proof.* Soundness: `fᵢ ≥ b′` on every segment, so `c ≥ b′` on the polygon, which lies in `Γ(y)`.
Completeness: choose `b″ ∈ (b′, E(y))` and a path `γ` with `min_γ c > b″` (definition of `sup`).
By uniform continuity pick points `γ(s₀=0), …, γ(s_m=1)` with consecutive distances `< η/2`,
`η := b″−b′`, and replace each interior point by a rational point within `η/4`, keeping `p₀ = y`
(if `y ∉ ℚ³` use the real certificate; the checker is stated for rational `y`) and `p_m` in the
bulk. Each point `q` of a segment `[p_k,p_{k+1}]` is within `η` of some `γ(s_k)`, so
`c(q) ≥ c(γ(s_k)) − η > b′` (F2a). ∎

## Theorem B3 (box-barrier certificates: soundness and completeness)
Grid: `g ∈ ℚ³`, `h > 0`, cells `C_k = g + h(k + [0,1]³)`. Let `A, S ⊂ ℤ³` be finite and disjoint with
(i) `y ∈ ⋃_{k∈A} C_k`; (ii) `k ∈ A ⇒ k ± e_d ∈ A ∪ S` for `d = 1,2,3`;
(iii) each `C_s` (`s ∈ S`) has all 8 corners in a single ball `B̄(x_{i(s)}, r_{i(s)} + b)`;
(iv) all corners `q` of all `C_k`, `k ∈ A`, satisfy `‖q−o‖ < R`.
*Soundness:* `E_X(y) ≤ b`. *Completeness:* if `c_X(y) > b` and `E_X(y) < b`, such a certificate
exists for every `h ≤ 2(b − E_X(y)) / (3√6)` (with bulk radius enlarged to `R + 2h√3`, allowed by
F1b), with `|A| + |S| = O(Leb(W ⊕ B(0,2h√3))/h³)` for the basin `W` of `y`.

*Proof of soundness.* Put `W_A = ⋃_{k∈A} C_k` (compact, ⊆ `B(o,R)` by (iv) and convexity of
cells). **Claim:** `∂W_A ⊆ ⋃_{s∈S} C_s`. Let `z ∈ ∂W_A` and suppose `z` lies in no shell cell.
The set of indices of cells containing `z` is `N(z) = N₁×N₂×N₃` with each `N_d` equal to
`{⌊(z_d−g_d)/h⌋}` or, when `(z_d−g_d)/h = m ∈ ℤ`, to `{m−1, m}`. Since `z ∈ W_A` (closed), some
`k ∈ A ∩ N(z)`. Any `k′ ∈ N(z)` is reached from `k` by changing one coordinate at a time by ±1
inside `N(z)`; each step goes from a cell of `A` to a face-neighbour containing `z`, which by (ii)
is in `A ∪ S` and, since it contains `z`, not in `S`; so `N(z) ⊆ A`. For `η` smaller than the
distance from `z` to every grid plane not containing `z`, the cube `{‖w−z‖_∞ < η}` lies in
`⋃_{k∈N(z)} C_k ⊆ W_A`, so `z ∈ int W_A`, contradicting `z ∈ ∂W_A`.
Now take `γ ∈ Γ(y)`. `γ(0) ∈ W_A` and `γ(1) ∉ W_A` (by (iv), `W_A ⊆ B(o,R)`). With
`s* = max{s : γ(s) ∈ W_A}` (closed set), `γ(s*) ∈ ∂W_A ⊆` some `C_s ⊆ B̄(x_{i(s)}, r_{i(s)}+b)`
(convexity: corners in a ball ⇒ cell in the ball). So `min_γ c ≤ f_{i(s)}(γ(s*)) ≤ b`. ∎

*Proof of completeness.* Let `e = E(y) < b`, `t = e + (b−e)/3`, and `W` the component of
`{c > t}` containing `y` (exists since `c(y) > b > t`); `W ∩ B_R = ∅` because `E(y) < t`. Let
`A = {k : C_k ∩ W ≠ ∅}` and `S = {k ∉ A : k is a face-neighbour of some element of A}`.
(i),(ii) hold by construction; (iv) holds with the enlarged bulk radius. For (iii), let `s ∈ S` be
adjacent to `k ∈ A`. `C_k ∪ C_s` is connected, meets `W`, and `C_s ∩ W = ∅`, so it contains a
point `z ∈ ∂W`, where `c(z) ≤ t`. Pick `i` with `fᵢ(z) = c(z)`. Every `q ∈ C_s` has
`‖q−z‖ ≤ diam(C_k ∪ C_s) = h√6`, so `fᵢ(q) ≤ t + h√6 ≤ b` when `h√6 ≤ 2(b−e)/3`. ∎

## Lemma B4a (power-covering lemma)
If closed balls `B̄(xᵢ,Rᵢ)`, `i ∈ t`, have a common point `p`, then
`conv{xᵢ : i ∈ t} ⊆ ⋃_{i∈t} B̄(xᵢ,Rᵢ)`.

*Proof.* For `q = Σ λᵢxᵢ` (`λ ∈ Δ`), expanding `‖p−xᵢ‖² = ‖p−q‖² + 2(p−q)·(q−xᵢ) + ‖q−xᵢ‖²` and
using `Σλᵢ(q−xᵢ) = 0` gives `Σλᵢ‖p−xᵢ‖² = ‖p−q‖² + Σλᵢ‖q−xᵢ‖²`. Hence
`Σλᵢ(‖q−xᵢ‖² − Rᵢ²) = Σλᵢ(‖p−xᵢ‖² − Rᵢ²) − ‖p−q‖² ≤ 0`, so some `i` with `λᵢ > 0` has
`‖q−xᵢ‖ ≤ Rᵢ`. ∎ (Formalized in Lean: `PocketGenesis/Covering.lean`.)

## Lemma B4b (mod-2 winding parity of a realized 2-cycle)
Let `T` be a finite set of 3-subsets of `[n]` forming a mod-2 cycle and `φ: [n] → ℝ³` a
realization; `|T|_φ = ⋃_{t∈T} conv φ(t)`. For `z ∉ |T|_φ` let `w_T(z) ∈ ℤ₂` be the mod-2 degree of
`F_z: |K_T| → S²`, `q ↦ (φ(q)−z)/‖φ(q)−z‖` (affine extension of `φ` on each triangle) evaluated on
the fundamental class `[T] ∈ H₂(|K_T|;ℤ₂)`. Then (a) for every direction `u` such that the ray
`z + ℝ_{>0}u` misses all realized edges and vertices and crosses every realized triangle it meets
transversally at a relative-interior point, `w_T(z)` equals the parity of the number of triangles
crossed; (b) `w_T` is constant along any path in `ℝ³ ∖ |T|_φ`; (c) `w_T(z) = 0` if `z` lies outside
`conv φ([n])` or on its boundary sphere-side (in particular for `z ∈ B_R`, `R ≥ ρ_X`); (d) `w_T` is
invariant under any homotopy `(φ_s, z_s)` with `z_s ∉ |T|_{φ_s}` for all `s`, and under isometries.

*Proof.* Standard mod-2 degree theory for piecewise-linear maps of a ℤ₂-pseudocycle (e.g.
Rourke–Sanderson, *Introduction to PL topology*; Matoušek, *Using the Borsuk–Ulam theorem*, §2.4):
(a) `u` is a regular value whose preimages are the transversal crossings, at each of which `F_z`
is a local homeomorphism; (b),(d) homotopy invariance of `F_*: H₂(|K_T|;ℤ₂) → H₂(S²;ℤ₂)`;
(c) a ray pointing away from the convex hull meets nothing. For (c) with `z` on the bulk sphere:
`|T|_φ ⊆ B̄(o,ρ_X)` and the outward radial ray from `z` leaves `B̄(o,ρ_X)` immediately. ∎

## Theorem B4 (Čech-cycle barrier certificates: soundness)
If `T` is a mod-2 cycle, each `t ∈ T` has a witness `p_t` with `‖p_t−xᵢ‖ ≤ rᵢ + b` (`i ∈ t`),
`c_X(y) > b`, and `w_T(y) = 1` (realization `φ(i) = xᵢ`), then `E_X(y) ≤ b`.

*Proof.* By B4a each realized triangle lies in `⋃_{i∈t} B̄(xᵢ,rᵢ+b) ⊆ {c ≤ b}`; so `y ∉ |T|`. If a
path `γ ∈ Γ(y)` avoided `|T|`, B4b(b,c) would give `1 = w_T(y) = w_T(γ(1)) = 0`. So `γ` meets
`|T| ⊆ {c ≤ b}` and `min_γ c ≤ b`. ∎

## Theorem B5 (Čech-cycle certificates are complete; rational witnesses exist)
Assume the weighted points `(xᵢ, (rᵢ+b)²)` are in general position (no degenerate power
configurations). If `c_X(y) > b` and `E_X(y) ≤ b`, there is a certificate as in B4 with `T` a set
of triangles of the weighted alpha complex `K_b` (the dual complex of `U_b = ⋃B̄(xᵢ,rᵢ+b)`); if
the data are rational, every `p_t` can be chosen rational.

*Proof.* (1) Edelsbrunner [*The union of balls and its dual shape*, DCG 13 (1995)]: `|K_b| ⊆ U_b`
and `|K_b|` is a deformation retract of `U_b`. (2) `E(y) ≤ b` and `c(y) > b` mean `y` lies in a
bounded component of `ℝ³ ∖ U_b` (F3). Alexander duality `H̃₀(S³∖A;ℤ₂) ≅ Ȟ²(A;ℤ₂)` for the
compact ENRs `A = |K_b| ⊆ U_b` is natural for inclusions; since `Ȟ²(U_b) → Ȟ²(|K_b|)` is an
isomorphism, so is `H̃₀(S³∖U_b) → H̃₀(S³∖|K_b|)`. Injectivity means `y` and `∞` lie in different
components of `S³ ∖ |K_b|`. (3) `|K_b|` is a subcomplex of the (embedded) weighted Delaunay
triangulation `D`; components of `ℝ³∖|K_b|` are unions of open simplices of `D ∖ K_b` (and the
outer region), and two tetrahedra of `D∖K_b` sharing a triangle not in `K_b` lie in the same
component. Let `Q` be the (bounded) component of `y` and `c_Q` the 3-chain of its tetrahedra. Then
`T := ∂c_Q` is a mod-2 cycle; each `σ ∈ T` has `Q` on exactly one side, hence `σ ∈ K_b` (a
triangle outside `K_b` would connect its two sides). For a point `z` interior to a tetrahedron of
`Q`, `w_T(z)` is the multiplicity of `z` in `c_Q`, i.e. 1; `y` has a neighbourhood in `Q` avoiding
`|T|`, so `w_T(y) = 1`. (4) `σ ∈ K_b` means `V_σ ∩ ⋂_{i∈σ} B̄(xᵢ,rᵢ+b) ≠ ∅` (`V_σ` the dual power
face). The minimum of the common power distance over `V_σ` is attained at the power orthocentre
of some coface `τ ⊇ σ` lying in `V_τ`; that point is rational (it solves linear equations with
rational coefficients) and lies in all balls of `σ`. ∎

*Remark.* Soundness (B4) needs **no** general position and accepts *any* Čech triangles (not only
alpha-complex ones), any self-intersecting realization, and any rational witnesses: the untrusted
solver may be floating-point; only the checker must be exact.

## Theorem R1 (robust certificates with fixed spatial witnesses)
Let `U = U(X⁰;ε,η)` and choose the bulk with `R ≥ maxᵢ(‖x⁰ᵢ−o‖ + εᵢ)`.
(a) `‖y−x⁰ᵢ‖ ≥ a + r⁰ᵢ + εᵢ + ηᵢ ∀i ⇒ c_X(y) ≥ a ∀X∈U`.
(b) `dist(x⁰ᵢ,[p_k,p_{k+1}]) ≥ b′ + r⁰ᵢ + εᵢ + ηᵢ ∀i,k` (and B2's endpoint condition)
`⇒ E_X(y) ≥ b′ ∀X ∈ U`.
(c) A box barrier whose shell cells lie in `B̄(x⁰_{i(s)}, r⁰_{i(s)} + b − ε_{i(s)} − η_{i(s)})`
`⇒ E_X(y) ≤ b ∀X ∈ U`.
*Proof.* Triangle inequalities: `‖z−xᵢ‖ − rᵢ ≥ ‖z−x⁰ᵢ‖ − εᵢ − r⁰ᵢ − ηᵢ` and
`‖z−xᵢ‖ − rᵢ ≤ ‖z−x⁰ᵢ‖ + εᵢ − r⁰ᵢ + ηᵢ`; then B2, B3 for each `X`. ∎

## Theorem R2 (robust moving-barrier certificates: one combinatorial cycle certifies a whole region)
Let `T` be a mod-2 cycle with witnesses `p_t` satisfying
`‖p_t − x⁰ᵢ‖ ≤ r⁰ᵢ + b − εᵢ − ηᵢ` (`i ∈ t`), let `a > b` and
`‖y − x⁰ᵢ‖ ≥ a + r⁰ᵢ + εᵢ + ηᵢ` for **all** `i`, and `w_T(y; X⁰) = 1`. Then every `X ∈ U(X⁰;ε,η)`
satisfies `c_X(y) ≥ a` and `E_X(y) ≤ b`; i.e. `P_{a,b}` holds on all of `U`, with the same witness.

*Proof.* For `X ∈ U` let `X(s) = X⁰ + s(X−X⁰)` (`U` is convex). For each `s`,
`‖p_t − xᵢ(s)‖ ≤ ‖p_t − x⁰ᵢ‖ + sεᵢ ≤ r⁰ᵢ − ηᵢ + b ≤ rᵢ(s) + b`, so every triangle of `T` is Čech at
level `b` for `X(s)` and, by B4a, `|T|_{X(s)} ⊆ {c_{X(s)} ≤ b}`; by R1(a), `c_{X(s)}(y) ≥ a > b`,
so `y ∉ |T|_{X(s)}`. By B4b(d), `w_T(y; X(s))` is constant, equal to 1. Apply B4 at `s = 1`. ∎

## Theorem R3 (intrinsic certificates: invariant under rigid motion)
Let `y(X) = Σ_k α_k x_k`, `p_t(X) = Σ_k β_{t,k} x_k` with `Σα = Σβ_t = 1`, and set
`ε_y = Σ|α_k|ε_k`, `ε_t = Σ_k |β_{t,k}| ε_k`. If `‖p_t(X⁰) − x⁰ᵢ‖ ≤ r⁰ᵢ + b − εᵢ − ηᵢ − ε_t`,
`‖y(X⁰) − x⁰ᵢ‖ ≥ a + r⁰ᵢ + εᵢ + ηᵢ + ε_y`, and `w_T(y(X⁰); X⁰) = 1`, then `P_{a,b}(gX)` holds with
witness `g·y(X) = y(gX)` for every `X ∈ U(X⁰;ε,η)` and every isometry `g`.
*Proof.* As R2 with the witnesses moving affinely along the homotopy (`‖y(X(s)) − y(X⁰)‖ ≤ ε_y`
etc.); isometry invariance of every condition, including the mod-2 parity, by B4b(d). ∎

*Consequence (certificate-based pocket identity).* An intrinsic certificate `(α, T, β)` is a
conformation-independent object; it can be checked on every frame of an ensemble without
alignment. We define two conformations to *share a certified pocket* if one intrinsic certificate
is valid on both. This gives a rigorous, alignment-free notion of pocket correspondence; Theorem
M1 shows that no *global* continuous labelling can exist, so local certificates are the right
granularity.

## Theorem R4 (monotone locality of confinement)
A Čech-cycle certificate `(T, {p_t}, y)` for `X` remains valid for every `X′` obtained from `X` by
adding, removing, moving or resizing atoms **not** in `V(T) = ⋃T`, provided every atom of `X′`
satisfies the chamber inequality `‖y − x′ᵢ‖ ≥ a + r′ᵢ`. Confinement is thus a property of the
cycle atoms alone; only the chamber clearance is a global (one `O(n)` scan) condition.
*Proof.* The witness inequalities and `w_T(y)` involve only atoms of `V(T)`; apply B4 to `X′`. ∎

## Theorem R5 (robust completeness with margin)
In general position, if `c_{X⁰}(y) ≥ a + ε̂` and `E_{X⁰}(y) ≤ b − ε̂ < c_{X⁰}(y)` with
`ε̂ = maxᵢ(εᵢ+ηᵢ)`, then certificates as in R2 and (for `E_{X⁰}(y) < b − ε̂`) R1(c) exist.
*Proof.* B5 (resp. B3) at level `b − ε̂` for `X⁰`; witnesses lie in balls of radius
`r⁰ᵢ + b − ε̂ ≤ r⁰ᵢ + b − εᵢ − ηᵢ`; the chamber inequality holds with margin `ε̂`. ∎

## Complexity of checking
With `n` atoms, a path certificate with `m` segments is checked in `O(nm)` exact operations (or
`O(m + n)` with a verified spatial index), a box barrier in `O(|A| + |S|)` plus `O(n)` for the
chamber, and a Čech-cycle certificate in `O(|T|)` for witnesses + `O(|T|)` ray crossings + `O(n)`
for the chamber. All arithmetic is in ℚ with bit-length linear in the input bit-length (squared
distances of rationals); no square roots are ever taken.
