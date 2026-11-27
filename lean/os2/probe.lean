import Mathlib
open MeasureTheory Topology
example : TopologicalSpace.PseudoMetrizableSpace ℝ := inferInstance
example : TopologicalSpace.PseudoMetrizableSpace (Set.Icc (-1 : ℝ) 1) := inferInstance
example (κ : Type*) : TopologicalSpace.PseudoMetrizableSpace (κ → Set.Icc (-1 : ℝ) 1) := inferInstance
example (κ : Type*) : CompactSpace (κ → Set.Icc (-1 : ℝ) 1) := inferInstance
example (κ : Type*) : T2Space (κ → Set.Icc (-1 : ℝ) 1) := inferInstance
