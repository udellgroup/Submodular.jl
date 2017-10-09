#############################################################################
# rank_of_graph_matroid.jl
# given a graph G = (V, E)
# f(S) = n - # of connected components of the subgraph of G induced by S
#############################################################################

export rank_of_graph_matroid
export sign, monotonicity, modularity, evaluate

type RankOfGraphicMatroidAtom <: CombiFunc
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractGraph, CombiSet}
  size::Tuple{Int, Int}
  graph::AbstractGraph
  setvariables::Array{CombiSet}

  function RankOfGraphicMatroidAtom(g::AbstractGraph, S::CombiSet)
    if nv(g) != S.cardinality
      error("Cannot define a modular function when the vector's size is different from the (card(baseset), 1).")
    else
      children = (g, S)
      setvariables = get_sv(S)
      return new(:rankofgraphmatroid, hash(children), children, (1, 1), g, setvariables)
    end
  end
end

rank_of_graph_matroid(g::AbstractGraph, S::CombiSet) = RankOfGraphicMatroidAtom(g, S)

function sign(x::RankOfGraphicMatroidAtom)
  return Positive()
end

function monotonicity(x::RankOfGraphicMatroidAtom)
  return Nonincreasing()
end

function modularity(x::RankOfGraphicMatroidAtom)
  return SubModularity()    # up to verification
end

function evaluate(x::RankOfGraphicMatroidAtom)
  sub_graph, vmap = induced_subgraph(x.children[1], get_elements(x.children[2]))
  return x.children[2].cardinality - length(connected_components(sub_graph))
end
