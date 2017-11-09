#############################################################################
# polyhedra.jl
# Handles constraints specifying the feasible region to be polyhedra
# associated with a submodular function f: 2ⱽ → ℝ, where V is the base set.
# The included polyhedra are:
#
## The submodular polyhedron: SubmodPoly = {s ∈ ℝᵖ, ∀ A ⊆ V, s(A) ≤ F(A)}.
## The base polyhedron: BasePoly = {s ∈ ℝᵖ, s(V) = F(V), ∀ A ⊆ V, s(A) ≤ F(A)}
##                               = SubmodPoly ∩ {s(V) = F(v)}.
## The positive polyhedron: PosPoly = {s ∈ ℝᵖ₊, ∀ A ⊆ V, s(A) ≤ F(A)}
##                                  = SubmodPoly ∩ ℝᵖ₊.
## The symmetric submodular polyhedron: SymPoly = {s ∈ ℝᵖ, ∀ A ⊆ V, |s|(A) ≤ F(A)}
##                                              = {s ∈ ℝᵖ, |s| ∈ SubmodPoly}.
# Here p = |V|, i.e. the cardinality of the base set.
#############################################################################

import Base.in

export AssocPoly
export SubmodPoly, BasePoly, PosPoly, SymPoly
export in, fenchel

abstract type AssocPoly <: ContiSet end

type SubmodPoly{T} <: AssocPoly
  head::Symbol
  F::T              # the submodular function
  V::AbstractArray # indexes of the base set
end

function SubmodPoly{T<:SubmodFunc}(F::T)
  @assert evaluate(F, [])[1] == 0
  x = get_sv(F)[1]
  V = x.baseset
  return(SubmodPoly(:subpoly, F, V))
end

type BasePoly{T} <: AssocPoly
  head::Symbol
  F::T
  V::AbstractArray
end

function BasePoly{T<:SubmodFunc}(F::T)
  @assert evaluate(F, [])[1] == 0
  x = get_sv(F)[1]
  V = x.baseset
  return(BasePoly(:basepoly, F, V))
end

type PosPoly{T} <: AssocPoly
  head::Symbol
  F::T
  V::AbstractArray
end

function PosPoly{T<:SubmodFunc}(F::T)
  @assert evaluate(F, [])[1] == 0
  x = get_sv(F)[1]
  V = x.baseset
  return(PosPoly(:pospoly, F, V))
end

type SymPoly{T} <: AssocPoly
  head::Symbol
  F::T
  V::AbstractArray
end

function SymPoly{T<:SubmodFunc}(F::T)
  @assert evaluate(F, [])[1] == 0
  x = get_sv(F)[1]
  V = x.baseset
  return(SymPoly(:sympoly, F, V))
end

function get_sv(p::AssocPoly)
  return get_sv(p.F)
end

function get_cv(p::AssocPoly)
  return Variable[]
end

# check if w is in p

function in(w::AbstractArray, p::AssocPoly, Tol::Float64 = 1e-3)
  n = length(w)
  @assert length(p.V) == n
  x = prox(p, w, Tol)
  if sum(abs.(x - w)) < Tol
    return true
  else
    return false
  end
end

# compute min w^x st x in p

# fenchel(p, w) = greedy(p.F, -w, p.V), w is nonpositve,
#               = -Inf,                 w is not nonpositive
function fenchel(p::SubmodPoly, w::AbstractArray)
  n = length(w)
  @assert length(p.V) == n
  if any(x -> x>0, w)
    return -Inf
  else
    return greedy(p.F, -w)
    # return greedy(p.F, -w, p.V)
  end
end

# fenchel(p, w) = greedy(p.F, -w, p.V)
function fenchel(p::BasePoly, w::AbstractArray)
  n = length(w)
  @assert length(p.V) == n
  return greedy(p.F, -w)
  # return greedy(p.F, -w, p.V)
end

# fenchel(p, w) = -sign(w_{-}) * greedy(p.F, -w_{-}, p.V)
function fenchel(p::PosPoly, w::AbstractArray)
  n = length(w)
  @assert length(p.V) == n
  u = 0.5 * (w - abs.(w))
  v = greedy(p.F, -u)
  # v = greedy(p.F, -u, p.V)
  v = -sign.(u) .* v
  return v
end

# fenchel(p, w) = -sign(w) .* greedy(p.F, |w|, p.V)
function fenchel(p::SymPoly, w::AbstractArray)
  n = length(w)
  @assert length(p.V) == n
  v = greedy(p.F, abs.(w))
  # v = greedy(p.F, abs.(w), p.V)
  v = -sign.(w) .* v
  return v
end
