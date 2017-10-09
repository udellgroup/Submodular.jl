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

function sign(x::CutAtom)
  return Positive()
end

function monotonicity(x::CutAtom)
  return NoSign()
end

function modularity(x::CutAtom)
  return SubModularity()    # up to verification
end

function evaluate(x::CutAtom)
  # weights = weights(x.children[1])
  cut = 0.0
  set = get_elements(x.children[2])
  for e in edges(x.children[1])
    if in(src(e), set) + in(dst(e), set) == 1
      cut += 1
    end
  end
  val = zeros(1, 1)
  val[1] = cut
  return val
end
