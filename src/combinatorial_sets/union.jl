#############################################################################
# union.jl
# Handles unions of combinatorial sets.
#############################################################################

import Base: union

export union, get_elements

type UnionAtom <: CombiSet
  head::Symbol
  id_hash::UInt64
  children::Tuple
  elements::ValOrNothing
  baseset::AbstractArray
  cardinality::Int
  value::ValOrNothing
  sign::Sign

  function UnionAtom(elements::AbstractArray; baseset = elements::AbstractArray,
  children = (elements, )::Tuple, sign = NoSign()::Sign)
    this = new(:union, 0, children, elements, baseset, length(baseset), nothing, NoSign())
    this.id_hash = object_id(this)
    return this
  end
end

function union(set::AllCombiSet...)
  sets = collect(set)
  n = length(sets)
  baseset = Set([])
  if n < 2
    error("Cannot compute the union of one set.")
  else
    elements = []
    if typeof(sets) <: Array{Number}
      for i = 1:n
        elements = vcat(elements, sets[i])
      end
      elements = collect(Set(elements))
      return elements
    else
      elements = Set([])
      for i = 1:n
        if typeof(sets[i]) <: CombiSet
          baseset = union(baseset, Set(sets[i].baseset))
          elements = union(elements, Set(get_elements(sets[i])))
        else
          seti = Set(sets[i])
          baseset = union(baseset, seti)
          elements = union(elements, seti)
        end
      end
    end
    baseset = collect(baseset)
    elements = collect(elements)
    newset = UnionAtom(elements, children = set, baseset = baseset, sign = Nondecreasing)
    return newset
  end
end

function get_elements(x::UnionAtom)
  elements = []
  n = length(x.children)
  for i = 1:n
    if typeof(x.children[i]) <: CombiSet
      elements = union(elements, Set(get_elements(x.children[i])))
    else
      seti = Set(x.children[i])
      elements = union(elements, seti)
    end
  end
  x.elements = collect(elements)
  return x.elements
end
