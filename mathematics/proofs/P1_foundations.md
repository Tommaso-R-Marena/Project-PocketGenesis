# Part I — Foundations of the clearance / escape-clearance model

Conventions as in `mathematics/definitions/DEFINITIONS.md`. `f_i(y) = ‖y−xᵢ‖−rᵢ`,
`c = minᵢ fᵢ`, bulk `B_R = {‖z−o‖ ≥ R}` with `R ≥ ρ_X = maxᵢ‖xᵢ−o‖`.

---

## Lemma F0 (clearance is distance to the union outside it)
For `y ∉ Ω(X)`, `c_X(y) = dist(y, Ω(X))`; for `y ∈ Ω(X)`, `c_X(y) ≤ 0`.

*Proof.* `dist(y, B̄(x,r)) = max(‖y−x‖−r, 0)`. If `y ∉ Ω` every `fᵢ(y) > 0`, so
`dist(y,Ω) = minᵢ max(fᵢ(y),0) = minᵢ fᵢ(y)`. If `y ∈ B̄(xᵢ,rᵢ)` then `fᵢ(y) ≤ 0`. ∎

## Lemma F1 (radial monotonicity; independence of the artificial boundary)
(a) If `‖z−o‖ ≥ ρ_X` and `u` is the unit vector `(z−o)/‖z−o‖` (any unit vector if `z=o`), then
`t ↦ c_X(z+tu)` is nondecreasing on `[0,∞)`.
(b) For `R′ ≥ R ≥ ρ_X`, the escape clearance computed with bulk `B_R` equals the one computed with
`B_{R′}`. In particular `E_X` does not depend on the artificial domain as long as the bulk sphere
encloses all atom **centres**.
(c) If `y ∈ B_R` then `E_X(y) = c_X(y)`.

*Proof.* (a) `‖z+tu−xᵢ‖² = ‖z−xᵢ‖² + 2t(z−xᵢ)·u + t²` and
`(z−xᵢ)·u = ‖z−o‖ − (xᵢ−o)·u ≥ ‖z−o‖ − ‖xᵢ−o‖ ≥ 0`; so each `fᵢ(z+tu)` is nondecreasing, hence
so is their minimum.
(b) `Γ_{R′}(y) ⊆ Γ_R(y)` gives `E^{(R)} ≥ E^{(R′)}`. Conversely, a path in `Γ_R(y)` ending at `z`
(`‖z−o‖ ≥ R ≥ ρ_X`) is extended by the radial segment from `z` to radius `R′`; by (a) the
minimum of `c` does not decrease, so `E^{(R′)} ≥ E^{(R)}`.
(c) The constant path is admissible, so `E(y) ≥ c(y)`; always `E ≤ c` (F2). ∎

*Remark (artificial boundaries).* (b) is sharp in the following sense: if `R < ρ_X` the bulk
may contain an atom centre and `E` can change (a path may "escape" by ending inside the protein).
Every computation in this project uses `R ≥ ρ_X`; certificates record `(o, R²)` and the checker
verifies `R² ≥ ρ_X²`.

## Lemma F2 (Lipschitz structure)
(a) `c_X` is 1-Lipschitz. (b) `E_X ≤ c_X`, hence `D_X ≥ 0`. (c) `E_X` is 1-Lipschitz on ℝ³.
(d) `D_X` is 2-Lipschitz.

*Proof.* (a) Each `fᵢ` is 1-Lipschitz; minima preserve this. (b) `γ(0) = y`. (c) Given `y, y′` and
any `γ ∈ Γ(y)`, the concatenation of the segment `[y′,y]` with `γ` lies in `Γ(y′)` and its minimum
is `≥ min(c(y) − ‖y−y′‖, min_γ c)`. Taking the supremum over `γ`,
`E(y′) ≥ min(c(y) − ‖y−y′‖, E(y)) ≥ E(y) − ‖y−y′‖` using `E(y) ≤ c(y)`. Swap roles.
(d) from (a), (c). ∎

## Lemma F3 (the supremum in E is attained)
For every `y` there is `γ* ∈ Γ(y)` with `min_s c_X(γ*(s)) = E_X(y)`. Consequently
`E_X(y) ≥ t ⇔ y` is joined to `B_R` by a path in `{c_X ≥ t}`, and
`E_X(y) > t ⇔ y` is joined to `B_R` by a path in `{c_X > t}`,
`E_X(y) ≤ t ⇔` no path from `y` to `B_R` stays in `{c_X > t}`.

*Proof.* Let `t* = E(y)`, `K_t = {c ≥ t} ∩ B̄(o,R)`, `S = ∂B(o,R)`. For each `k` there is a path
with minimum `≥ t* − 1/k`; truncating it at its first hitting time of `S` gives a path in
`K_{t*−1/k}` from `y` to `S`. Let `C_k` be the connected component of `K_{t*−1/k}` containing `y`;
then `C_k ∩ S ≠ ∅`, `C_{k+1} ⊆ C_k`, and each `C_k` is a continuum (compact connected). The
intersection `C = ⋂C_k` is a continuum containing `y`, meets `S` (nested nonempty compacts
`C_k ∩ S`), and lies in `K_{t*}`. `K_{t*}` is semialgebraic (`fᵢ ≥ t ⇔ (t+rᵢ ≤ 0) ∨
(‖·−xᵢ‖² ≥ (t+rᵢ)²)`), and connected components of semialgebraic sets are semialgebraically path
connected [Bochnak–Coste–Roy, *Real Algebraic Geometry*, Prop. 2.5.13 & Thm. 2.4.5]. So the
component of `K_{t*}` containing `C` contains a path from `y` to a point of `S ⊆ B_R`, along which
`c ≥ t*`. The equivalences follow: the third is the negation of the second, and the second uses
only the definition of `sup` (no attainment needed). ∎

## Lemma F4 (monotonicity and invariances)
(a) If `X ⊆ X′` (more atoms) or `X′` is `X` with radii increased, then `c_{X′} ≤ c_X` and
`E_{X′} ≤ E_X` (common bulk). (b) For any isometry `g` of ℝ³, `c_{gX}(gy) = c_X(y)`,
`E_{gX}(gy) = E_X(y)`. (c) `c`, `E`, `D`, `P_{a,b}` depend only on the multiset of balls.

*Proof.* (a) pointwise smaller integrand over the same path family. (b) `g` maps `Γ_X(y)` onto
`Γ_{gX}(gy)` (with bulk `gB_R`, admissible by F1(b)). (c) definitions. ∎

## Lemma F5 (margin characterisation)
For `a > b`: `τ_{a,b}(X) ≥ 0 ⇔ P_{a,b}(X)`, and the supremum defining `τ` is a maximum over
`B̄(o,R)`.

*Proof.* On `B_R`, `E = c` (F1c) so `min(c−a, b−E) ≤ (b−a)/2 < 0` there. On the compact ball the
function `y ↦ min(c(y)−a, b−E(y))` is continuous (F2) and attains its maximum `m`. If `m ≥ 0` the
maximiser is a chamber witness; if `m < 0` then `τ = max(m, sup over B_R) < 0` and no witness
exists. ∎

## Theorem F6 (basin decomposition of the trapped region)
Let `W` be a connected component of the open set `{D_X > 0}`. Then
1. `E_X` is constant on `W`; call the value `e_W`;
2. `W` is a bounded connected component of `{c_X > e_W}`, and `W ∩ B_R = ∅`;
3. for every `y` with `D_X(y) > 0`, the basin containing `y` is the component of
   `{c_X > E_X(y)}` containing `y`;
4. (gate) `∂W` contains a point `s` with `c_X(s) = e_W` that is joined to `B_R` by a path in
   `{c_X ≥ e_W}`;
5. `E_X = c_X` on the complement of `⋃ W`.

*Proof.* (1) *Local constancy on `{E < c}`.* If `E(y) < c(y)`, put `δ = (c(y)−E(y))/2` and let
`‖y′−y‖ < δ`. On `[y,y′]`, `c > E(y)`. Concatenating, `E(y′) ≥ min(min_{[y,y′]} c, E(y)) = E(y)`;
and if `E(y′) > E(y)` the same concatenation from `y` would give `E(y) ≥ min(min c, E(y′)) > E(y)`.
So `E(y′) = E(y)`. A locally constant function is constant on connected sets.
(2) `W ⊆ {c > e_W}` and `W` is connected, so `W ⊆ C`, the component of `{c > e_W}` containing it.
If `C` met `B_R`, for `y ∈ W` there would be a path in `C` (open connected ⇒ path connected) from
`y` to `B_R` with (compactness) minimum `> e_W`, contradicting `E(y) = e_W`. Hence `C ⊆ B(o,R)`.
For `y′ ∈ C`: joining `y′` to `y ∈ W` inside `C` and then following an optimal escape path of `y`
(F3) shows `E(y′) ≥ e_W`; if `E(y′) > e_W`, the reverse concatenation would give `E(y) > e_W`.
So `E ≡ e_W` on `C`, `D > 0` on `C`, `C ⊆ {D > 0}` connected, hence `C = W`.
(3) is (2) with `e_W = E(y)`. (4) Let `γ*` be an optimal escape path of `y ∈ W` (F3). The set
`{s : γ*(s) ∉ W}` is closed and nonempty; its minimum `s₁ > 0` gives `γ*(s₁) ∈ ∂W`; since
`∂W ∩ {c > e_W} = ∅` (a boundary point of a component of an open set is not in the set) and
`c ≥ e_W` on `γ*`, we get `c(γ*(s₁)) = e_W`, and `γ*|_{[s₁,1]}` stays in `{c ≥ e_W}`.
(5) `D = 0` there. ∎

*Interpretation.* Pockets "at all probe scales simultaneously" are canonically the basins: the
bounded components of superlevel sets of the clearance that merge **directly into the bulk**
(children of the root in the merge tree). A sub-chamber nested inside a larger confined region
belongs to the larger basin; the full hierarchy is the merge tree (Part IV). The trapping
functional `D` therefore measures *confinement relative to the bulk*, not concavity: a hemispherical
bowl has `D ≡ 0` (Counterexample CE-F1 in `COUNTEREXAMPLES.md`).

## Proposition F7 (what D distinguishes — and what it does not)
In the union-of-balls model:
(i) a wide shallow depression has small `D` (chamber and mouth clearances nearly equal);
(ii) a deep narrow-necked pocket has `D ≥ c(chamber) − c(neck)` large;
(iii) a sealed cavity has `E < 0`, so `D > c ≥ 0`, distinguishing it from open pockets by the sign
of `E` (not by `D`);
(iv) a tunnel with several exits: `E` sees only the widest exit; the number of exits is invisible
to `E` and `D` and is recorded instead by `H₁` classes (Part IV, D3.5);
(v) a highly concave but unconfined region has `D = 0` though its buriedness is high.
*Proof.* (i)-(iii),(v) are immediate from F6 and the definitions; (iv): two configurations
differing only by closing a narrower second exit have identical `E` (F3: the widest exit realises
the supremum) but different `H₁({c > b})` for `b` below the second exit's bottleneck. ∎
