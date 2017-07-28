#############################################################################
# polyhedra.jl
# Handles constraints specifying the feasible region to be polyhedra
# associated with a submodular function f: 2ⱽ → ℝ, where V is the base set.
# The included polyhedra are:
## The submodular polyhedron: SubPoly = {s ∈ ℝᵖ, ∀ A ⊆ V, s(A) ≤ F(A)}.
## The base polyhedron: BasePoly = {s ∈ ℝᵖ, s(V) = F(V), ∀ A ⊆ V, s(A) ≤ F(A)}
##                               = subpoly ∩ {s(V) = F(v)}.
## The positive polyhedron: PosPoly = {s ∈ ℝᵖ₊, ∀ A ⊆ V, s(A) ≤ F(A)}
##                                  = subpoly ∩ ℝᵖ₊.
## The symmetric submodular polyhedron: SymPoly = {s ∈ ℝᵖ, ∀ A ⊆ V, |s|(A) ≤ F(A)}
##                                              = {s ∈ ℝᵖ, |s| ∈ subpoly}.
# Here p = |V|, i.e. the cardinality of the base set.
#############################################################################

import Base.in
export SubPoly, BasePoly, PosPoly, SymPoly
export in, fenchel

type SubPoly
  f::Function # the submodular function
  V::AbstractVector # indexes of the base set

  function SubPoly(f::Function, V::AbstractVector)
    @assert f([]) == 0
    new(f, V)
  end

  function SubPoly(f::Function, n::Int)
    @assert f([]) == 0
    new(f, collect(1:n))
  end
end

type BasePoly # throw errors if f is not defined on V
  f::Function
  V::AbstractVector

  function BasePoly(f::Function, V::AbstractVector)
    @assert f([]) == 0
    new(f, V)
  end

  function BasePoly(f::Function, n::Int)
    @assert f([]) == 0
    new(f, collect(1:n))
  end
end

type PosPoly # throw errors if f is not defined on V
  f::Function
  V::Array{Int}

  function PosPoly(f::Function, V::AbstractVector)
    @assert f([]) == 0
    new(f, V)
  end

  function PosPoly(f::Function, n::Int)
    @assert f([]) == 0
    new(f, collect(1:n))
  end
end

type SymPoly # throw errors if f is not defined on V
  f::Function
  V::AbstractVector

  function SymPoly(f::Function, V::AbstractVector)
    @assert f([]) == 0
    new(f, V)
  end

  function SymPoly(f::Function, n::Int)
    @assert f([]) == 0
    new(f, collect(1:n))
  end
end

# check if u is in p

function in(w::AbstractVector, p::SubPoly)
  n = length(w)
  @assert length(p.V) == n
  checker = true
  ordering = sortperm(w, rev = true)
  h = 0
  ind = []
  for i = 1:length(w)
    h += - p.f(p.V[ind]) - w[i] + p.f(p.V[push!(ind, ordering[i])])
    if h < 0
      checker = false
      break
    end
  end
  return checker
end

function in(w::AbstractVector, p::BasePoly)
  n = length(w)
  @assert length(p.V) == n
  checker = true
  ordering = sortperm(w, rev = true)
  h = 0
  ind = []
  for i = 1:length(w)
    h += - p.f(p.V[ind]) - w[i] + p.f(p.V[push!(ind, ordering[i])])
    if h < 0
      checker = false
      break
    end
  end
  if checker == true
    checker = (sum(w) == p.f(p.V))
  end
  return checker
end

function in(w::AbstractVector, p::PosPoly)
  n = length(w)
  @assert length(p.V) == n
  checker = true
  ordering = sortperm(w, rev = true)
  h = 0
  ind = []
  for i = 1:length(w)
    h += - p.f(p.V[ind]) - w[i] + p.f(p.V[push!(ind, ordering[i])])
    if h < 0
      checker = false
      break
    end
  end
  if checker == true
    checker = all(x -> x>=0, w)
  end
  return checker
end

function in(u::AbstractVector, p::SymPoly)
  w = abs.(u)
  n = length(w)
  @assert length(p.V) == n
  checker = true
  ordering = sortperm(w, rev = true)
  h = 0
  ind = []
  for i = 1:length(w)
    h += - p.f(p.V[ind]) - w[i] + p.f(p.V[push!(ind, ordering[i])])
    if h < 0
      checker = false
      break
    end
  end
  return checker
end

# compute min c^x st p
# fenchel(p, u) = -lovaszext(p.f, p.V, -u), u is nonpositve,
#               = -Inf,                     u is not nonpositive
function fenchel(p::SubPoly, w::AbstractVector)
  n = length(w)
  @assert length(p.V) == n
  if any(x -> x>0, w)
    return -Inf
  else
    return -lovaszext(p.f, -w, p.V)
  end
end

# fenchel(p, u) = -lovaszext(p.f, p.V, -u)
function fenchel(p::BasePoly, w::AbstractVector)
  n = length(w)
  @assert length(p.V) == n
  return -lovaszext(p.f, -w, p.V)
end

# fenchel(p, u) = -lovaszext(p.f, p.V, -u_{+})
function fenchel(p::PosPoly, w::AbstractVector)
  n = length(w)
  @assert length(p.V) == n
  u = 0.5 * (w - abs.(w))
  return -lovaszext(p.f, -u, p.V)
end

# fenchel(p, u) = -lovaszext(p.f, p.V, -|u|)
function fenchel(p::SymPoly, w::AbstractVector)
  n = length(w)
  @assert length(p.V) == n
  return -lovaszext(p.f, abs.(w), p.V)
end
