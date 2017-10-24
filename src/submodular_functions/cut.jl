#############################################################################
# cut.jl
# given a graph G = (V, E)
# f(S) = \sum of the weights of the edges spanning across S and (V\S)
#############################################################################

export cut
export sign, monotonicity, modularity, evaluate

type CutAtom <: CombiFunc
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractGraph, CombiSet}
  size::Tuple{Int, Int}
  graph::AbstractGraph
  setvariables::Array{CombiSet}

  function CutAtom(g::AbstractGraph, S::CombiSet)
    if nv(g) != S.cardinality
      error("Cannot define a modular function when the vector's size is different from the (card(baseset), 1).")
    else
      children = (g, S)
      setvariables = get_sv(S)
      return new(:cut, hash(children), children, (1, 1), g, setvariables)
    end
  end
end

cut(g::AbstractGraph, S::CombiSet) = CutAtom(g, S)

function sign(F::CutAtom)
  return Positive()
end

function monotonicity(F::CutAtom)
  return (NoMonotonicity(), )
end

function modularity(F::CutAtom)
  w = weights(F.children[1])
  if all(x -> x>=0, w)
    return SubModularity()
  elseif all(x -> x <=0, w)
    return SuperModularity()
  end
end

function evaluate(F::CutAtom)
  w = weights(F.children[1])
  cut = 0.0
  set = get_elements(F.children[2])
  for e in edges(F.children[1])
    u, v = src(e), dst(e)
    if in(u, set) + in(v, set) == 1
      cut += w[u, v]
    end
  end
  return cut
end
