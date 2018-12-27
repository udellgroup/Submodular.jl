#############################################################################
# union.jl
# Handles unions of combinatorial sets.
#############################################################################

import Base: union

export union, get_elements

mutable struct UnionAtom <: CombiSet
  head::Symbol
  id_hash::UInt64
  children::Tuple
  elements::ValOrNothing
  baseset::AbstractArray
  cardinality::Int

  function UnionAtom(elements::ValOrNothing, children::Tuple, baseset::AbstractArray)
    this = new(:union, 0, children, elements, baseset, length(baseset))
    this.id_hash = object_id(this)
    return this
  end
end

function union(set::AllCombiSet...)
  sets = collect(set)
  baseset = Set([])
  n = length(sets)
  if n < 2
    error("Cannot compute the union of one set.")
  else
    if isa(sets, Array{Number})
      elements = []
      for i = 1:n
        elements = vcat(elements, sets[i])
      end
      elements = collect(Set(elements))
      return elements
    else
      for i = 1:n
        if isa(sets[i], CombiSet)
          baseset = union(baseset, Set(sets[i].baseset))
        else
          seti = Set(sets[i])
          baseset = union(baseset, seti)
        end
      end
      baseset = collect(baseset)
      newset = UnionAtom(nothing, set, baseset)
      return newset
    end
  end
end

function get_elements(x::UnionAtom)
  elements = []
  n = length(x.children)
  for i = 1:n
    if isa(x.children[i], CombiSet)
      elements = union(elements, Set(get_elements(x.children[i])))
    else
      seti = Set(x.children[i])
      elements = union(elements, seti)
    end
  end
  x.elements = collect(elements)
  return x.elements
end
