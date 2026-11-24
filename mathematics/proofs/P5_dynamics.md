# Part V — Dynamic pocket formation: certified events, finiteness, monodromy

Trajectory `X: [0,T] → Conf_n`, per-atom speeds `vᵢ` (D5.1), `v = maxᵢ vᵢ`, common bulk with
`R ≥ sup_t ρ_{X(t)}`. Write `τ(t) = τ_{a,b}(X(t))`.

## Theorem Dy1 (Lipschitz margins in time)
`t ↦ c_{X(t)}(y)`, `t ↦ E_{X(t)}(y)` and `t ↦ τ(t)` are `v`-Lipschitz; `t ↦ D_{X(t)}(y)` is
`2v`-Lipschitz. With per-atom speeds, the local versions of S2/S3 hold with `εᵢ = vᵢ|t−t′|`.
*Proof.* S1 and S2/S3 applied to `X(t), X(t′)`. ∎

## Theorem Dy2 (certified event localisation from finitely many samples)
Let `t_k < t_{k+1}`, `Δ = t_{k+1} − t_k`, and certified bounds `τ(t_k) ∈ [ℓ_k, u_k]`.
(a) *No event:* if `u_k + u_{k+1} + vΔ < 0` then `τ < 0` (pocket absent) on `[t_k, t_{k+1}]`; if
`ℓ_k + ℓ_{k+1} − vΔ > 0` then `τ > 0` (pocket present) throughout.
(b) *Certified formation:* if `u_k < 0 < ℓ_{k+1}`, the pocket set `{t : τ(t) ≥ 0}` meets
`[t_k,t_{k+1}]`, its first point `t_f` satisfies `t_k + |u_k|/v ≤ t_f ≤ t_{k+1} − ℓ_{k+1}/v`, and the
last absence time before `t_f` lies in the same interval. Symmetrically for closure.
(c) Otherwise the interval is reported `UNRESOLVED`; refinement (an extra sample) is needed.
*Proof.* `τ(t) ≤ min(u_k + v(t−t_k), u_{k+1} + v(t_{k+1}−t))`, whose maximum over the interval is
`(u_k + u_{k+1} + vΔ)/2`. (b) Intermediate value theorem for the continuous `τ`, and the Lipschitz
bound forbids `τ ≥ 0` within `|u_k|/v` of `t_k` and forces `τ > 0` within `ℓ_{k+1}/v` of
`t_{k+1}`. ∎ (Lean: `PocketGenesis/Dynamics.lean`, theorem `no_event_of_margins` and
`event_interval`.)

## Proposition Dy3 (sampling complexity of event-free certification)
In the oracle model (the algorithm sees only exact values `τ(t_k)` and the bound `v`):
(a) a sample set (all `τ(t_k) < 0`) certifies `τ < 0` on `[0,T]` iff the open intervals
`(t_k − |τ(t_k)|/v, t_k + |τ(t_k)|/v)` cover `[0,T]`: at an uncovered `t*` the fooling function
`g(t) = min_k (τ(t_k) + v|t−t_k|)` — the largest `v`-Lipschitz function agreeing with the samples —
satisfies `g(t*) ≥ 0`, so a configuration path realising `g` cannot be excluded;
(b) the greedy rule `t_{k+1} = t_k + |τ(t_k)|/v` certifies with at most
`1 + 2∫₀ᵀ v/|τ(t)| dt` samples; (c) every certifying set has at least `vT / (2 max|τ|)` samples.
*Proof.* (a) Sufficiency: `τ(t) ≤ τ(t_k) + v|t−t_k| < 0` on the `k`-th interval. Necessity: `g`
is `v`-Lipschitz, `g(t_j) = τ(t_j)` (as `τ(t_k) + v|t_j−t_k| ≥ τ(t_j)`), and at an uncovered `t*`,
`v|t*−t_k| ≥ |τ(t_k)|` for all `k`, so `g(t*) ≥ 0`. (The statement is about the oracle model;
whether every `v`-Lipschitz margin function is realised by a trajectory is not claimed.)
(b) On `[t_k,t_{k+1}]`, `|τ| ≤ 2|τ(t_k)|`, so `∫ v/|τ| ≥ ½` per step. (c) covering length. ∎
(Status: standard Lipschitz certification, adapted; see `NOVELTY_AUDIT.md`.)

## Theorem Dy4 (finitely many events for tame motions; infinitely many for smooth ones)
(a) If `X(t)` is semialgebraic in `t` (e.g. piecewise polynomial interpolation), then `τ` is a
semialgebraic function on `[0,T]`; the pocket set `{t : P_{a,b}(X(t))}` is a finite union of points
and intervals, there are finitely many events, and every event time is a real algebraic number.
The same finiteness holds for globally subanalytic trajectories.
(b) There is a `C^∞` trajectory (with bounded speed) with infinitely many events in `[0,1]`.
*Proof.* (a) `c` is semialgebraic in `(t,y)`. The relation "`y` and `∂B(o,R)` lie in the same
connected component of `K_{t,s} = {z ∈ B̄(o,R) : c_{X(t)}(z) ≥ s}`" is semialgebraic in `(t,s,y)`
(connected components of a semialgebraic family are uniformly semialgebraic — Hardt's triviality
theorem, BCR Thm 9.3.2), so by F3 `E_{X(t)}(y) = max{s : …}` is semialgebraic, and so is `τ`
(a sup over `y` of a semialgebraic function). A one-variable semialgebraic function has finitely
many sign changes, and its zeros are algebraic. The subanalytic case uses the o-minimality of
`ℝ_an` identically. (b) Take the octahedral cage (Example S1 without `q`) with radii
`r(t) = r₀ + e^{−1/t²} sin(1/t)` (`r(0) = r₀`), `r₀ = ρ√(2/3) − b`, and `a < ρ − r₀ − 1`.
Then `E_{X(t)}(0) = ρ√(2/3) − r(t)` crosses `b` at every zero of `sin(1/t)`. ∎

*Consequence (minimum temporal information).* Without a speed bound, no finite set of samples
certifies the absence of events (Dy3a with `v = ∞`); with a speed bound, certification is
possible iff `τ` is bounded away from zero, at a cost governed by `∫ v/|τ|`. Tameness alone
(e.g. "the motion is polynomial of unknown degree") gives finiteness but no a-priori bound.

## Theorem M1 (pocket monodromy: no global continuous pocket identity)
There is a smooth loop `L: S¹ → Conf_n` (a *peristaltic torus*) and scales `a > b`, along which
`m = 3` pockets persist with no events, such that continuation around the loop permutes them
cyclically. Consequently there is no labelling of pockets by a fixed finite set that is
(i) locally constant under event-free continuation and (ii) bijective at every configuration of
the loop; in particular, no global, continuous pocket identity exists on any subset of
`Conf_n` containing `L`.

*Construction.* A torus tube: `N` rings of `k` atoms (radius `r`), ring `i` centred on the core
circle of radius `R_c` at angle `θᵢ = 2πi/N`, its atoms on a circle of radius
`ρᵢ(φ) = ρ₀ − A·g(mθᵢ − φ)` in the ring's normal plane, where `g` is a smooth `2π`-periodic bump
(`g = 1` near `0`, `g = 0` outside a small arc). Constricted rings (`ρᵢ` small) plug the lumen.
Since `ρᵢ(φ+2π) = ρᵢ(φ)`, `L` is a closed loop in *labelled* configuration space; the `m`
constrictions sit at `θ ≡ (φ + 2πj)/m`, so after one period each constriction — and each chamber
between consecutive constrictions — has advanced by `2π/m`.

*Proof.* (1) Symmetry: `L(φ + 2πm/N)` is the rotation by `2π/N` about the torus axis of `L(φ)`
with ring labels shifted by one, so it suffices to certify the pocket structure for
`φ ∈ [0, 2πm/N]`. (2) On a partition of this interval into sub-intervals `I`, robust intrinsic
certificates (R2/R3, per-atom tolerance `εᵢ = sup|ρ̇ᵢ|·|I|/2`) certify for every `φ ∈ I`: three
chamber witnesses `y₁,y₂,y₃` with `c ≥ a`, three Čech cycles `T_j` with `w_{T_j}(y_j) = 1` and
`w_{T_j}(y_l) = 0` for `l ≠ j` (so the three chambers lie in distinct basins), and path
certificates joining the witnesses of adjacent sub-intervals inside `{c > b}` (continuation).
These are checked exactly by the trusted checker (`experiments/synthetic/monodromy_torus.py`;
certificate files in `certificates/monodromy/`). (3) Following the certified witnesses around the
loop, the chamber starting between constrictions `j` and `j+1` ends between `j+1` and `j+2`.
(4) If a labelling with (i),(ii) existed, the pocket covering over `S¹` would be trivial, but its
monodromy is a 3-cycle (transitive). ∎

*Status.* Computer-assisted proof: the mathematical reduction (1),(3),(4) is on paper; step (2) is
an exact rational computation verified by the trusted checker. The physical reading — peristaltic
transporters and rotary mechanisms (functional rotation in multidrug efflux pumps) move pockets
along non-contractible loops of conformation space — is an interpretation, not a claim about
any specific protein.
