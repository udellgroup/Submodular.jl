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
  f::T              # the submodular function
  V::AbstractArray # indexes of the base set
end

function SubmodPoly{T<:CombiFunc}(f::T)
  @assert evaluate(f, [])[1] == 0
  x = get_sv(f)[1]
  V = x.baseset
  return(SubmodPoly(:subpoly, f, V))
end

type BasePoly{T} <: AssocPoly
  head::Symbol
  f::T
  V::AbstractArray
end

function BasePoly{T<:CombiFunc}(f::T)
  @assert evaluate(f, [])[1] == 0
  x = get_sv(f)[1]
  V = x.baseset
  return(BasePoly(:basepoly, f, V))
end

type PosPoly{T} <: AssocPoly
  head::Symbol
  f::T
  V::AbstractArray
end

function PosPoly{T<:CombiFunc}(f::T)
  @assert evaluate(f, [])[1] == 0
  x = get_sv(f)[1]
  V = x.baseset
  return(PosPoly(:pospoly, f, V))
end

type SymPoly{T} <: AssocPoly
  head::Symbol
  f::T
  V::AbstractArray
end

function SymPoly{T<:CombiFunc}(f::T)
  @assert evaluate(f, [])[1] == 0
  x = get_sv(f)[1]
  V = x.baseset
  return(SymPoly(:sympoly, f, V))
end

function get_sv(p::AssocPoly)
  return get_sv(p.f)
end

# check if w is in p

function in(p::AssocPoly, w::AbstractArray)
  n = length(w)
  @assert length(p.V) == n
  x = prox(p, w)
  if sum(abs.(x - w)) < TOL
    return true
  else
    return false
  end
end

# compute min w^x st x in p

# fenchel(p, w) = greedy(p.f, -w, p.V), w is nonpositve,
#               = -Inf,                 w is not nonpositive
function fenchel(p::SubmodPoly, w::AbstractArray)
  n = length(w)
  @assert length(p.V) == n
  if any(x -> x>0, w)
    return -Inf
  else
    return greedy(p.f, -w)
    # return greedy(p.f, -w, p.V)
  end
end

# fenchel(p, w) = greedy(p.f, -w, p.V)
function fenchel(p::BasePoly, w::AbstractArray)
  n = length(w)
  @assert length(p.V) == n
  return greedy(p.f, -w)
  # return greedy(p.f, -w, p.V)
end

# fenchel(p, w) = -sign(w_{-}) * greedy(p.f, -w_{-}, p.V)
function fenchel(p::PosPoly, w::AbstractArray)
  n = length(w)
  @assert length(p.V) == n
  u = 0.5 * (w - abs.(w))
  v = greedy(p.f, -u)
  # v = greedy(p.f, -u, p.V)
  v = -sign.(u) .* v
  return v
end

# fenchel(p, w) = -sign(w) .* greedy(p.f, |w|, p.V)
function fenchel(p::SymPoly, w::AbstractArray)
  n = length(w)
  @assert length(p.V) == n
  v = greedy(p.f, abs.(w))
  # v = greedy(p.f, abs.(w), p.V)
  v = -sign.(w) .* v
  return v
end
