#############################################################################
# cut.jl
# given a graph G = (V, E)
# f(S) = \sum of the edge_weights of the edges spanning across S and (V\S)
#############################################################################

export cut
export sign, monotonicity, modularity, evaluate

mutable struct CutAtom <: SubmodFunc
  head::Symbol
  id_hash::UInt64
  children::Tuple{WeightedGraph, CombiSet}
  size::Tuple{Int, Int}
  graph::WeightedGraph
  setvariables::Array{CombiSet}

  function CutAtom(g::WeightedGraph, S::CombiSet)
    if nv(g) != S.cardinality
      error("Cannot define a cut function when the number of vertices in the graph is different from the size of the set variable.")
    else
      children = (g, S)
      setvariables = get_sv(S)
      return new(:cut, hash(children), children, (1, 1), g, setvariables)
    end
  end
end

cut(g::WeightedGraph, S::CombiSet) = CutAtom(g, S)

function sign(F::CutAtom)
  return NoSign()
end

function monotonicity(F::CutAtom)
  return (NoMonotonicity(), )
end

function modularity(F::CutAtom)
  w = values(edge_weights(F.children[1]))
  if all(x -> x>=0, w)
    return SubModularity()
  elseif all(x -> x <=0, w)
    return SuperModularity()
  elseif all(x -> x==0, w)
    return ConstModularity()
  else
    return NotDetermined()
  end
end

function evaluate(F::CutAtom)
  w = edge_weights(F.children[1])
  cut = 0.0
  set = get_elements(F.children[2])
  for e in edges(F.children[1])
    u, v = src(e), dst(e)
    if in(u, set) + in(v, set) == 1
      cut += w[(u, v)]
    end
  end
  return cut
end
