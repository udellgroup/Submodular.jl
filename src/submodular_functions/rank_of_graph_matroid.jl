#############################################################################
# rank_of_graph_matroid.jl
# given a graph G = (V, E)
# f(S) = n - # of connected components of the subgraph of G induced by S
#############################################################################

export rank_of_graph_matroid
export sign, monotonicity, modularity, evaluate

type RankOfGraphicMatroidAtom <: SubmodFunc
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractGraph, CombiSet}
  size::Tuple{Int, Int}
  graph::AbstractGraph
  setvariables::Array{CombiSet}

  function RankOfGraphicMatroidAtom(g::AbstractGraph, S::CombiSet)
    if ne(g) != S.cardinality
      error("Cannot define a rank function of a graphic matroid when the number of edges in the graph is different from the size of the set variable.")
    else
      children = (g, S)
      setvariables = get_sv(S)
      return new(:rankofgraphmatroid, hash(children), children, (1, 1), g, setvariables)
    end
  end
end

rank_of_graph_matroid(g::AbstractGraph, S::CombiSet) = RankOfGraphicMatroidAtom(g, S)

function sign(F::RankOfGraphicMatroidAtom)
  return Positive()
end

function monotonicity(F::RankOfGraphicMatroidAtom)
  return (Nonincreasing(), )
end

function modularity(F::RankOfGraphicMatroidAtom)
  return SubModularity()    # up to verification
end

function evaluate(F::RankOfGraphicMatroidAtom)
  elist = AbstractEdge[]
  i = 1
  set = get_elements(F.children[2])
  for e in edges(F.children[1])
    if in(i, set) == true
      u, v = src(e), dst(e)
      push!(elist, Edge(u, v))
    end
    i += 1
  end
  sub_graph, vmap = induced_subgraph(F.children[1], elist)
  val = nv(sub_graph) - length(connected_components(sub_graph))
  return val
end
