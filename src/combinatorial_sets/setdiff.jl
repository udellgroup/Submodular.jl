#############################################################################
# setdiff.jl
# Handles set differences of combinatorial sets.
#############################################################################

import Base: setdiff

export setdiff, get_elements

type SetDiffAtom <: CombiSet
  head::Symbol
  id_hash::UInt64
  children::Tuple
  elements::ValOrNothing
  baseset::AbstractArray
  cardinality::Int
  value::ValOrNothing
  sign::Sign

  function SetDiffAtom(elements::AbstractArray; baseset = elements::AbstractArray,
  children = (elements, )::Tuple, sign = NoSign()::Sign)
    this = new(:setdiff, 0, children, elements, baseset, length(baseset), nothing, NoSign())
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
  set01 = Set(get_elements(set1))
  set02 = Set(set2)
  elements = collect(Set(setdiff(set01, set02)))
  newset = SetDiffAtom(elements, children = (set1, Constant(set2)), baseset = set1.baseset)
  return newset
end

function setdiff(set1::AbstractArray, set2::CombiSet)
  set01 = Set(set1)
  set02 = Set(get_elements(set2))
  elements = collect(Set(setdiff(set01, set02)))
  newset = SetDiffAtom(elements, children = (Constant(set1), set2), baseset = set2.baseset)
  return newset
end

function setdiff(set1::CombiSet, set2::CombiSet)
  set01 = Set(get_elements(set1))
  set02 = Set(get_elements(set2))
  elements = collect(Set(setdiff(set01, set02)))
  baseset = collect(union(Set(set1.baseset), Set(set2.baseset)))
  newset = SetDiffAtom(elements, children = (set1, set2), baseset = baseset)
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
