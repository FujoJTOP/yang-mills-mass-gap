import Mathlib
open MeasureTheory Topology
example (κ : Type*) [Countable κ] : TopologicalSpace.PseudoMetrizableSpace (κ → Set.Icc (-1 : ℝ) 1) := inferInstance
example (κ : Type*) [Countable κ] : SecondCountableTopology (κ → Set.Icc (-1 : ℝ) 1) := inferInstance
example (κ : Type*) [Countable κ] : BorelSpace (κ → Set.Icc (-1 : ℝ) 1) := inferInstance
example (κ : Type*) [Countable κ] : CompactSpace (ProbabilityMeasure (κ → Set.Icc (-1 : ℝ) 1)) := inferInstance
