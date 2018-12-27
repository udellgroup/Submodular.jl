#############################################################################
# setdiff.jl
# Handles set differences of combinatorial sets.
#############################################################################

import Base: setdiff

export setdiff, get_elements

mutable struct SetDiffAtom <: CombiSet
  head::Symbol
  id_hash::UInt64
  children::Tuple
  elements::ValOrNothing
  baseset::AbstractArray
  cardinality::Int

  function SetDiffAtom(elements::ValOrNothing, children::Tuple, baseset::AbstractArray)
    this = new(:setdiff, 0, children, elements, baseset, length(baseset))
    this.id_hash = object_id(this)
    return this
  end
end

function setdiff(set1, set2::AbstractArray)
  set01 = Set(set1)
  set02 = Set(set2)
  return collect(setdiff(set01, set02))
end

function setdiff(set1::CombiSet, set2::AbstractArray)
  newset = SetDiffAtom(nothing, (set1, Constant(set2)), set1.baseset)
  return newset
end

function setdiff(set1::AbstractArray, set2::CombiSet)
  newset = SetDiffAtom(nothing, (Constant(set1), set2), set2.baseset)
  return newset
end

function setdiff(set1::CombiSet, set2::CombiSet)
  baseset = collect(union(Set(set1.baseset), Set(set2.baseset)))
  newset = SetDiffAtom(nothing, (set1, set2), baseset)
  return newset
end

function get_elements(x::SetDiffAtom)
  if isa(x.children[1], Constant)
    set1 = Set(x.children[1].value)
    set2 = x.children[2].elements
  elseif isa(x.children[2], Constant)
    set1 = x.children[1].elements
    set2 = Set(x.children[2].value)
  else
    set1 = x.children[1].elements
    set2 = x.children[2].elements
  end
  x.elements = setdiff(set1, set2)
  return x.elements
end
