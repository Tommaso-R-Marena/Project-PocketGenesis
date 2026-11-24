# Part IV — The pocket diagram (relative superlevel persistence) and channels

Status: this part *organizes* the theory using established persistence machinery
(Edelsbrunner–Letscher–Zomorodian 2002; Cohen-Steiner–Edelsbrunner–Harer 2007; Chazal et al.
2009). The identification of pocket predicates with quadrant counts is elementary; we do not
claim it as a novel theorem (see `NOVELTY_AUDIT.md`).

## Theorem H1 (pocket predicates are rank invariants of a relative persistence module)
Let `S_t = {c_X ≥ t}` and `V_t = H₀(S_t ∪ B_R, B_R; ℤ₂)` with the maps `V_t → V_{t′}` (`t ≥ t′`)
induced by inclusion.
(a) `dim V_t` is the number of components of `S_t` not meeting `B_R` ("confined chambers at
level t"); `V` is q-tame (indeed of finite type) because `c_X` is semialgebraic.
(b) For `a > b`, `rank(V_a → V_b)` is the number of components of `S_b ∪ B_R` not containing
`B_R` that meet `S_a`; in particular
`rank(V_a → V_b) ≥ 1 ⇔ ∃y : c_X(y) ≥ a ∧ E_X(y) < b`.
(c) (Stability) `d_B(Dgm V^X, Dgm V^Y) ≤ ‖c_X − c_Y‖_∞ ≤ d_∞(X,Y)`.
(d) (Basins ↔ bars) Under genericity (distinct critical values), the bars whose death is a merge
into the bulk component are exactly `{(m_W, e_W) : W basin}`; all other bars are sub-chambers
dying into older confined components.

*Proof.* (a) Relative `H₀` of a pair counts components of the space not meeting the subspace
(`B_R` is connected). Semialgebraic sets have finitely many components, and the family `S_t` has
finitely many topological types (Hardt), giving finite type. (b) The image of `V_a → V_b` is
spanned by classes of components of `S_b ∪ B_R` that contain a component of `S_a ∪ B_R` not
meeting `B_R`; for `y` with `c(y) ≥ a`, `y`'s component in `S_b ∪ B_R` avoids `B_R` iff no path
in `{c ≥ b}` joins `y` to `B_R` iff `E(y) < b` (F3). (c) `S_t(X) ⊆ S_{t−ε}(Y)` and the bulk is
common, so the modules are `ε`-interleaved; apply the algebraic stability theorem. (d) Elder rule:
within a basin `W` the component born at the highest maximum `m_W` is the oldest, so it survives
every internal merge and dies exactly when `W` joins the bulk at `e_W` (F6). ∎

*Convention note.* The closed superlevel module gives the strict mouth inequality `E < b`; the
open version `{c > t}` gives `E ≤ b`. The two differ only at critical levels.

## Proposition H2 (channels are H₁ classes; mouths beyond the first are H₁)
At probe level `b`, `H₁({c > b}; ℤ₂) ≅ H₁(U_b; ℤ₂)` (Alexander duality in `S³`), each class
corresponding to a handle of the inflated protein (a ring of atoms around a tunnel). For a basin
with `k` mouths whose bottleneck values exceed `b`, the union of the basin with the bulk contains
`k − 1` independent loops (enter through one mouth, leave through another) that are nonzero in
`H₁({c > b})`. The persistence of such a class in the superlevel filtration is the interval
between the ring-closure value (an index-1 critical value of `c`) and the channel bottleneck (an
index-2 critical value).
*Proof sketch.* Alexander duality for compact polyhedral `U_b ⊂ S³`; the loop through two mouths
links the ring of atoms separating them. The critical-point indices follow from the Morse theory
of distance functions (Gershkovich–Rubinstein; Lieutier 2004; flow complex of Giesen–John). ∎
(Status: standard; recorded for completeness of the taxonomy, not as a contribution.)
