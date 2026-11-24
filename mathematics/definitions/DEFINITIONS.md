# Foundational definitions

All statements in this project are about the explicit mathematical models defined here.
Physical interpretation is discussed separately in `SCIENTIFIC_LIMITATIONS.md`.

Notation: `‖·‖` is the Euclidean norm on ℝ³; `B̄(x,ρ)` the closed ball, `B(x,ρ)` the open ball;
`S²` the unit sphere; `Leb` Lebesgue measure.

---

## D1. Atomistic configurations

**D1.1 (Configuration).** A *configuration* is a finite labelled family
`X = ((x₁,r₁),…,(xₙ,rₙ))`, `xᵢ ∈ ℝ³`, `rᵢ ≥ 0`. The label set `[n]` is part of the data; two
configurations of the *same molecule* share the label set (atom correspondence is given, never
inferred). `Conf_n = (ℝ³ × ℝ_{≥0})ⁿ`.

**D1.2 (Displacement norms).** For `X, Y ∈ Conf_n` let `εᵢ(X,Y) = ‖xᵢ − yᵢ‖ + |rᵢ − r′ᵢ|`
(per-atom displacement) and `d_∞(X,Y) = maxᵢ εᵢ(X,Y)`. A *per-atom uncertainty region* around
`X⁰` with tolerances `ε = (εᵢ)`, `η = (ηᵢ)` is
`U(X⁰; ε, η) = { X : ‖xᵢ − x⁰ᵢ‖ ≤ εᵢ, |rᵢ − r⁰ᵢ| ≤ ηᵢ, rᵢ ≥ 0 ∀i }` (a compact convex set).

**D1.3 (Occupied region).** `Ω(X) = ⋃ᵢ B̄(xᵢ, rᵢ)` (union-of-balls model; closed).

*Modeling choices made explicit:* radii are effective (e.g. van der Waals or united-atom) and are
inputs; hydrogens may be absent; coordinates carry experimental error. The union-of-balls model is
**not** the solvent-excluded surface, nor a density or energy model; see D7 for the generalization
to arbitrary potential fields, under which every escape-functional result of Part I–III holds.

## D2. Clearance, escape clearance, trapping depth

**D2.1 (Clearance field).** `f_i(y) = ‖y − xᵢ‖ − rᵢ`, `c_X(y) = minᵢ f_i(y)`.
For `y ∉ Ω(X)`, `c_X(y) = dist(y, Ω(X)) > 0`; for `y ∈ Ω(X)`, `c_X(y) ≤ 0`. A probe ball of radius
`ρ ≥ 0` centred at `y` is *disjoint from the open atoms* iff `c_X(y) ≥ ρ` and *disjoint from the
closed atoms* iff `c_X(y) > ρ`. Contact (`= ρ`) is treated explicitly in every statement.

**D2.2 (Probe-accessible region).** `F_ρ(X) = { y : c_X(y) > ρ }` (open), and its closed variant
`F̄_ρ(X) = { y : c_X(y) ≥ ρ }`.

**D2.3 (Bulk).** Fix a reference point `o` and a radius `R ≥ ρ_X := maxᵢ ‖xᵢ − o‖` (only centres
matter). The bulk is `B_R = { z : ‖z − o‖ ≥ R }`. Lemma F1 shows that every quantity below is
independent of `(o, R)` subject to `R ≥ ρ_X`; when comparing configurations we use a common `o`
and `R ≥ max(ρ_X, ρ_Y)`.

**D2.4 (Escape paths).** `Γ(y) = { γ ∈ C([0,1], ℝ³) : γ(0) = y, γ(1) ∈ B_R }`.

**D2.5 (Escape clearance).** `E_X(y) = sup_{γ ∈ Γ(y)} min_{s ∈ [0,1]} c_X(γ(s))`.
(The inner minimum exists by compactness; the supremum is attained, Lemma F3.) Interpretation: the
largest probe radius that can be moved from `y` to bulk without overlapping atoms; the optimal
value is the *bottleneck* (maximin, "widest path") value.

**D2.6 (Trapping depth).** `D_X(y) = c_X(y) − E_X(y) ≥ 0`.

**D2.7 (Margins).** For scales `a > b` the *pocket margin* is
`τ_{a,b}(X) = sup_y min( c_X(y) − a , b − E_X(y) )`.

## D3. Pocket predicates and taxonomy

Three scales are used, each explicit: a **chamber scale** `a` (radius of the largest probe the
chamber must hold), a **mouth scale** `b < a` (the mouth must not pass a probe larger than `b`),
and a **solvent scale** `w ≥ 0` (e.g. water, 1.4 Å).

**D3.1 (Confined pocket).** `P_{a,b}(X)` holds iff `∃ y : c_X(y) ≥ a ∧ E_X(y) ≤ b`.
Any such `y` is a *chamber witness*. Equivalently (Lemma F5) `τ_{a,b}(X) ≥ 0`.

**D3.2 (Pocket basins).** The *basins* of `X` are the connected components of `{ D_X > 0 }`.
Theorem F6: each basin `W` is a bounded connected component of `{ c_X > e_W }` on which
`E_X ≡ e_W` (the *mouth value*); `m_W := sup_W c_X` is the *chamber value*; `m_W − e_W` is the
*basin depth*. Basins are the children of the bulk in the merge tree of superlevel sets of `c_X`.

**D3.3 (Taxonomy at scales (a, b, w)).** Let `y` be a point with `c_X(y) ≥ a`.
- *Sealed cavity:* `E_X(y) < 0` (no probe of any radius ≥ 0 escapes; the wall has no crack).
- *Buried cavity (solvent-inaccessible):* `E_X(y) ≤ w` (a solvent probe cannot reach `y`).
- *Open confined pocket:* `w < E_X(y) ≤ b` (solvent reaches `y`, a ligand-size probe of radius
  `> b` cannot pass the mouth, and the chamber holds radius `a`).
- *Deep pocket of depth ≥ τ:* `D_X(y) ≥ τ`.
- *Unconfined region:* `E_X(y) = c_X(y)` (`D = 0`): a probe at `y` of the local maximal size can
  leave. Concavity is **not** confinement: a hemispherical bowl is unconfined (`D ≡ 0`); see
  *buriedness* (D3.6) for the complementary, non-topological notion.
- *Surface depression:* unconfined (`D_X(y) < τ_min`) but with buriedness `β_{b,L}(y) ≥ θ`.

**D3.4 (Mouth, gate, bottleneck).** For a basin `W` with mouth value `e_W`, a *gate* is a point
`s ∈ ∂W` with `c_X(s) = e_W` through which `W` touches the closure of the unbounded component of
`{c_X > e_W}`. Generically (Part VIII) a gate is an index-2 critical point of `c_X` whose *active
atoms* `I(s) = { i : f_i(s) = c_X(s) }` (three, generically) are the *gate atoms*; the triangle of
their centres is the *mouth window*. For a channel/tunnel the *bottleneck* is the maximin value of
`c_X` along it.

**D3.5 (Channels).** At probe scale `b`, a *through-channel* is a nonzero class of
`H₁({c_X > b}; ℤ₂)` (a loop of accessible space that cannot be contracted in accessible space).
By Alexander duality it corresponds to a handle of the inflated union `U_b = {c_X ≤ b}`. A pocket
with `k` mouths open at scale `b` contributes `k − 1` independent channel classes.

**D3.6 (Directional buriedness).** For `L > 0`,
`β_{b,L}(y) = σ({ u ∈ S² : [y, y+Lu] ∩ U_b ≠ ∅ }) / 4π` where `σ` is surface measure.
Necessary condition for confinement (Theorem C2): `E_X(y) ≤ b ⇒ β_{b,∞}(y) = 1`.

**D3.7 (Interfacial pocket).** For a complex `X = A ∪ B` (disjoint label sets), an *interfacial
chamber witness* is `y` with `c_X(y) ≥ a`, `E_X(y) ≤ b`, but `E_A(y) > b` and `E_B(y) > b`
(neither partner alone confines `y`; by monotonicity, Lemma F4, `c_A, c_B ≥ c_X ≥ a`).

**D3.8 (Transient and cryptic pockets).** Let `M ⊆ Conf_n` be an admissible family and `X⁰ ∈ M`.
- *Transient* (along a trajectory `X(t)`): the set `{ t : P_{a,b}(X(t)) }` is neither empty nor
  everything.
- *Cryptic relative to `(X⁰, M)`*: `¬P_{a,b}(X⁰)` and `∃ X ∈ M : P_{a,b}(X)`.
  Six **distinct** propositions are kept apart (they need different evidence):
  (G) geometrically possible: `M = Conf_n`; (S) sterically feasible: `M` = configurations
  satisfying explicit non-bonded exclusion and bonded constraints; (E) energetically admissible:
  `M = { U ≤ U_max }`; (T) thermodynamically populated: `π(A) ≥ p_min` for the canonical measure
  `π` (Part VII); (P) pathway-reachable: `X` is in the path component of `X⁰` in `M`;
  (K) kinetically relevant: reached within a time scale under a stated dynamics.
  This project proves results about (G), (S) in explicit models, (E)/(T) in explicit energy
  models, (P) as a complexity statement, and makes **no** claims about (K).

**D3.9 (Minimal formation displacement and work).**
`δ*(X⁰; M, a, b) = inf { d(X⁰, X) : X ∈ M, P_{a,b}(X) }` for a stated metric `d`
(`d_∞`, weighted per-atom norms, or Mahalanobis); and, for a harmonic model with stiffness `H`,
`W*(X⁰) = inf { ½ ξᵀHξ : P_{a,b}(X⁰ + ξ) }` (minimal elastic work).

## D4. Certificates (checkable witnesses)

All certificate data are rational (`ℚ`); all checks are exact.

**D4.1 (Chamber certificate).** A point `y ∈ ℚ³` and level `a ∈ ℚ` with
`‖y − xᵢ‖² ≥ (a + rᵢ)²` and `a + rᵢ ≥ 0` for all `i` (certifies `c_X(y) ≥ a`).

**D4.2 (Path certificate, lower bound on E).** A polygon `y = p₀, …, p_m ∈ ℚ³` with
`‖p_m − o‖² ≥ R²` and `dist(xᵢ, [p_k, p_{k+1}])² ≥ (b′ + rᵢ)²` for all `i, k` (with `b′ + rᵢ ≥ 0`).
Certifies `E_X(y) ≥ b′`.

**D4.3 (Box-barrier certificate, upper bound on E).** A grid `(g, h)`, `g ∈ ℚ³`, `h ∈ ℚ_{>0}`,
cells `C_k = g + h(k + [0,1]³)`, `k ∈ ℤ³`; finite disjoint sets `A, S ⊂ ℤ³` such that
(i) `y ∈ ⋃_{k∈A} C_k`; (ii) every face-neighbour `k ± e_d` of every `k ∈ A` lies in `A ∪ S`;
(iii) every `C_s`, `s ∈ S`, has all 8 corners in some single `B̄(xᵢ, rᵢ + b)`;
(iv) all corners of all `C_k`, `k ∈ A`, satisfy `‖· − o‖² < R²`. Certifies `E_X(y) ≤ b`.

**D4.4 (Čech-cycle barrier certificate, upper bound on E).** A finite set `T` of 3-subsets of
`[n]` that is a mod-2 cycle (every 2-subset lies in an even number of members of `T`); for each
`t ∈ T` a witness `p_t ∈ ℚ³` with `‖p_t − xᵢ‖² ≤ (rᵢ + b)²` for `i ∈ t`; a rational ray direction
`u` from `y` in general position w.r.t. the realized triangles, crossing an odd number of them; and
`c_X(y) > b`. Certifies `E_X(y) ≤ b`.

**D4.5 (Robust and intrinsic certificates).** Versions of D4.1–D4.4 valid for every `X` in an
uncertainty region `U(X⁰; ε, η)` (radii shrunk/inflated by the tolerances; Theorem R1, R2), and
versions where `y` and the witnesses are affine combinations of atom positions (invariant under
rigid motions; Theorem R3).

**D4.6 (Verdicts).** A checker returns `PROVED` (a valid certificate for the claim), `REFUTED`
(a valid certificate for the negation), or `UNKNOWN`. `UNKNOWN` is never converted.

## D5. Dynamics

**D5.1** A *trajectory* is a continuous `X: [0,T] → Conf_n`. It is *`v`-Lipschitz* if
`εᵢ(X(t), X(t′)) ≤ vᵢ |t − t′|` for all `i`. An *event* of `P_{a,b}` is a time at which the truth
value of `P_{a,b}(X(t))` is not locally constant.

**D5.2** A *pocket covering* over a loop `L: S¹ → Conf_n` without events is the set of pairs
`(φ, W)` with `W` a basin of `L(φ)` containing a chamber witness; *monodromy* is the permutation
of basins obtained by continuation around the loop (Theorem M1).

## D6. Thermodynamics

Collective-coordinate model: compact box `Θ ⊂ ℝᵏ`, continuous `X: Θ → Conf_n`, continuous
`U: Θ → ℝ`, reference measure `μ` (Lebesgue × continuous density), canonical
`π(dθ) = Z⁻¹ e^{−βU(θ)} μ(dθ)`. The pocket set `A = { θ : P_{a,b}(X(θ)) }` is closed (Lemma T0).
Harmonic model: `X = X⁰ + ξ`, `ξ ~ N(0, Σ)`, `Σ = k_B T · H⁺` (rigid-body modes removed).

## D7. Generalized escape functionals

For any continuous `V: ℝ³ → ℝ` (e.g. a probe–protein interaction energy), define the *entry
barrier* `B_V(y) = inf_{γ ∈ Γ(y)} max_s V(γ(s))`. Hard-sphere geometry is the special case
`V = −c_X` (`B_V = −E_X`). Duality, stability (`‖B_V − B_W‖_∞ ≤ ‖V − W‖_∞`), basin decomposition
and the path/box-barrier certificates hold verbatim for any continuous `V` (with "inside a ball"
replaced by any exact sufficient test for `V ≥ level` on a cell). This is the precise sense in
which geometric and energetic accessibility are two instances of one theory; it does **not** make
geometric accessibility imply energetic accessibility.
