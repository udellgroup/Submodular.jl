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

  function IntersectAtom(elements::ValOrNothing, children::Tuple, baseset::AbstractArray)
    this = new(:intersect, 0, children, elements, baseset, length(baseset))
    this.id_hash = object_id(this)
    return this
  end
end

function intersect(set::AllCombiSet...)
  baseset = Set([])
  sets = collect(set)
  n = length(sets)
  if n < 2
    error("Cannot compute the union of one set.")
  else
    if isa(sets, Array{Number})
      elements = Sets[i]
      for i = 1:n
        elements = intersect(elements, Set(sets[i]))
      end
      elements = collect(Set(elements))
      return elements
    else
      for i = 2:n
        if isa(sets[i], CombiSet)
          baseset = union(baseset, Set(sets[i].baseset))
        else
          seti = Set(sets[i])
          baseset = union(baseset, seti)
        end
      end
      baseset = collect(baseset)
      newset = IntersectAtom(nothing, set, baseset)
      return newset
    end
  end
end

function get_elements(x::IntersectAtom)
  n = length(x.children)
  if isa(x.children[1], CombiSet)
    elements = Set(x.children[1].elements)
  else
    elements = x.children[1]
  end
  for i = 2:n
    if isa(x.children[i], CombiSet)
      elements = intersect(elements, Set(get_elements(x.children[i])))
    else
      seti = Set(x.children[i])
      elements = intersect(elements, seti)
    end
  end
  x.elements = collect(elements)
  return x.elements
end
