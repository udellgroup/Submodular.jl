#############################################################################
# intersect.jl
# Handles intersections of combinatorial sets.
#############################################################################

import Base: intersect

export intersect, get_elements

type IntersectAtom <: CombiSet
  head::Symbol
  id_hash::UInt64
  children::Tuple
  elements::ValOrNothing
  baseset::AbstractArray
  cardinality::Int
  value::ValOrNothing
  sign::Sign

  function IntersectAtom(elements::AbstractArray; baseset = elements::AbstractArray,
  children = (elements, )::Tuple, sign = NoSign()::Sign)
    this = new(:intersect, 0, children, elements, baseset, length(baseset), nothing, NoSign())
    this.id_hash = object_id(this)
    return this
  end
end

function intersect(set::AllCombiSet...)
  sets = collect(set)
  baseset = Set([])
  sets = collect(set)
  n = length(sets)
  if n < 2
    error("Cannot compute the union of one set.")
  else
    if typeof(sets) <: Array{Number}
      elements = Sets[i]
      for i = 1:n
        elements = intersect(elements, Set(sets[i]))
      end
      elements = collect(Set(elements))
      return elements
    else
      if typeof(sets[1]) <: CombiSet
        elements = Set(sets[1].elements)
      else
        elements = sets[1]
      end
      for i = 2:n
        if typeof(sets[i]) <: CombiSet
          baseset = union(baseset, Set(sets[i].baseset))
          elements = intersect(elements, Set(get_elements(sets[i])))
        else
          seti = Set(sets[i])
          baseset = union(baseset, seti)
          elements = intersect(elements, seti)
        end
      end
    end
    baseset = collect(baseset)
    elements = collect(elements)
    newset = IntersectAtom(elements, children = set, baseset = baseset, sign = Nondecreasing)
    return newset
  end
end

function get_elements(x::IntersectAtom)
  n = length(x.children)
  if typeof(x.children[1]) <: CombiSet
    elements = Set(x.children[1].elements)
  else
    elements = x.children[1]
  end
  for i = 2:n
    if typeof(x.children[i]) <: CombiSet
      elements = intersect(elements, Set(get_elements(x.children[i])))
    else
      seti = Set(x.children[i])
      elements = intersect(elements, seti)
    end
  end
  x.elements = collect(elements)
  return x.elements
end
