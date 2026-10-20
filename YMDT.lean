import Mathlib

open Matrix
open scoped Matrix BigOperators

attribute [local instance] Matrix.linftyOpNormedAddCommGroup Matrix.linftyOpNormedSpace

namespace YMDT

/-- Smoothness order C^∞ as a `WithTop ℕ∞` value, avoiding the `∞` notation scope. -/
abbrev Sm : WithTop ℕ∞ := ((⊤ : ℕ∞) : WithTop ℕ∞)

/-! ### 1. Compact Lie group G ≤ U(N) and Lie algebra 𝔤 -/

abbrev Mat (N : ℕ) := Matrix (Fin N) (Fin N) ℂ

structure CompactMatrixLieGroup (N : ℕ) where
  carrier  : Subgroup (unitaryGroup (Fin N) ℂ)
  isClosed : IsClosed (carrier : Set (unitaryGroup (Fin N) ℂ))

namespace CompactMatrixLieGroup
variable {N : ℕ}

def Mem (G : CompactMatrixLieGroup N) (g : Mat N) : Prop :=
  ∃ u : unitaryGroup (Fin N) ℂ, u ∈ G.carrier ∧ (u : Mat N) = g

def lieAlg (G : CompactMatrixLieGroup N) : Set (Mat N) :=
  {X | ∃ c : ℝ → Mat N, c 0 = 1 ∧ (∀ t, G.Mem (c t)) ∧ HasDerivAt c X 0}

def bracket (X Y : Mat N) : Mat N := X * Y - Y * X

end CompactMatrixLieGroup

/-! ### 2. Torus 𝕋⁴ with periodic lift -/

abbrev E4 := Fin 4 → ℝ
abbrev Tor4 := Fin 4 → AddCircle (1 : ℝ)

def toTor (x : E4) : Tor4 := fun i => ((x i : ℝ) : AddCircle (1 : ℝ))

def Periodic4 {α : Type*} (f : E4 → α) : Prop :=
  ∀ (ν : Fin 4) (x : E4), f (x + Pi.single ν (1 : ℝ)) = f x

noncomputable section

def pd {N : ℕ} (μ : Fin 4) (f : E4 → Mat N) (x : E4) : Mat N :=
  deriv (fun s : ℝ => f (x + s • Pi.single μ (1 : ℝ))) 0

/-! ### 3. Connections 𝒜 and gauge group 𝒢 -/

abbrev RawConn (N : ℕ) := E4 → Fin 4 → Mat N
abbrev RawGauge (N : ℕ) := E4 → Mat N

structure IsConn {N : ℕ} (G : CompactMatrixLieGroup N) (A : RawConn N) : Prop where
  smooth    : ∀ μ : Fin 4, ContDiff ℝ Sm (fun x => A x μ)
  periodic  : Periodic4 A
  valuesLie : ∀ x μ, A x μ ∈ G.lieAlg

structure IsGauge {N : ℕ} (G : CompactMatrixLieGroup N) (g : RawGauge N) : Prop where
  smooth   : ContDiff ℝ Sm g
  periodic : Periodic4 g
  valuesG  : ∀ x, G.Mem (g x)

def Conn {N : ℕ} (G : CompactMatrixLieGroup N) := {A : RawConn N // IsConn G A}
def CalG {N : ℕ} (G : CompactMatrixLieGroup N) := {g : RawGauge N // IsGauge G g}

/-! ### 4. Gauge action A^g = g⁻¹ A g + g⁻¹ dg -/

def gaugeAct {N : ℕ} (A : RawConn N) (g : RawGauge N) : RawConn N :=
  fun x μ => (g x)ᴴ * A x μ * g x + (g x)ᴴ * pd μ g x

def DirectGaugeRel {N : ℕ} (G : CompactMatrixLieGroup N) (A B : Conn G) : Prop :=
  ∃ g : CalG G, B.1 = gaugeAct A.1 g.1

abbrev ConnQuot {N : ℕ} (G : CompactMatrixLieGroup N) :=
  Quot (Relation.EqvGen (DirectGaugeRel G))

/-! ### 5. Curvature, covariant derivative, d*, drift ψ = -d*B -/

def covD {N : ℕ} (A : RawConn N) (μ : Fin 4) (φ : RawGauge N) (x : E4) : Mat N :=
  pd μ φ x + A x μ * φ x - φ x * A x μ

def curv {N : ℕ} (A : RawConn N) (μ ν : Fin 4) (x : E4) : Mat N :=
  pd μ (fun y => A y ν) x - pd ν (fun y => A y μ) x + A x μ * A x ν - A x ν * A x μ

def dStarF {N : ℕ} (A : RawConn N) (μ : Fin 4) (x : E4) : Mat N :=
  -∑ ν : Fin 4, covD A ν (fun y => curv A ν μ y) x

def dStar {N : ℕ} (B : RawConn N) (x : E4) : Mat N :=
  -∑ μ : Fin 4, pd μ (fun y => B y μ) x

def psi {N : ℕ} (B : RawConn N) (x : E4) : Mat N := -dStar B x

/-! ### 6. Parallel transport, Wilson loop trace -/

structure IsClosedLoop (γ : ℝ → E4) : Prop where
  smooth : ContDiff ℝ Sm γ
  period : ∃ n : Fin 4 → ℤ, γ 1 = γ 0 + (fun i => ((n i : ℤ) : ℝ))

def IsHol {N : ℕ} (A : RawConn N) (γ : ℝ → E4) (U : ℝ → Mat N) : Prop :=
  U 0 = 1 ∧ ∀ t, HasDerivAt U
    (-(∑ μ : Fin 4, (deriv (fun s => γ s μ) t : ℝ) • A (γ t) μ) * U t) t

def WilsonTr {N : ℕ} (_A : RawConn N) (_γ : ℝ → E4) (U : ℝ → Mat N) : ℂ :=
  (U 1).trace

def IsYMFlow {N : ℕ} (A : ℝ → RawConn N) : Prop :=
  ∀ t x μ, HasDerivAt (fun s => A s x μ) (-dStarF (A t) μ x) t

def IsDeTurckFlow {N : ℕ} (B : ℝ → RawConn N) : Prop :=
  ∀ t x μ, HasDerivAt (fun s => B s x μ)
    (-dStarF (B t) μ x + covD (B t) μ (psi (B t)) x) t

def IsDriftGauge {N : ℕ} (A : ℝ → RawConn N) (h : ℝ → RawGauge N) : Prop :=
  ∀ t x, HasDerivAt (fun s => h s x)
    (h t x * psi (gaugeAct (A t) (h t)) x) t

/-! ### Lemma 2.2 -/

def Lemma22 : Prop :=
  (∀ {N : ℕ} (G : CompactMatrixLieGroup N)
    (B : RawConn N), IsConn G B → ∀ (h : RawGauge N), IsGauge G h →
    ∀ (γ : ℝ → E4), IsClosedLoop γ → ∀ (U V : ℝ → Mat N),
      IsHol B γ U → IsHol (gaugeAct B h) γ V →
        WilsonTr B γ U = WilsonTr (gaugeAct B h) γ V)
  ∧
  (∀ {N : ℕ} (G : CompactMatrixLieGroup N)
    (A : ℝ → RawConn N) (h : ℝ → RawGauge N),
      (∀ t, IsConn G (A t)) → (∀ t, IsGauge G (h t)) →
      IsYMFlow A → IsDriftGauge A h →
      IsDeTurckFlow (fun t => gaugeAct (A t) (h t)) ∧
      ∀ (t : ℝ) (γ : ℝ → E4), IsClosedLoop γ → ∀ U V : ℝ → Mat N,
        IsHol (A t) γ U → IsHol (gaugeAct (A t) (h t)) γ V →
          WilsonTr (A t) γ U = WilsonTr (gaugeAct (A t) (h t)) γ V)

end

end YMDT
