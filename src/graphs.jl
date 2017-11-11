#############################################################################
# weights_of_graphs.jl
# handles the weights of AbstractGraph, complementary of LightGraphs.
#############################################################################

export WeightedGraph
export add_edge!, weights
export nv, ne, edges, induced_subgraph, connected_components

type WeightedGraph
  weights::Dict{Tuple{Int64,Int64},Number}
  graph::AbstractGraph

  function WeightedGraph(graph::AbstractGraph)
    return new(Dict{Tuple{Int64,Int64}, Number}(), graph)
  end
end

WeightedGraph(n::Int) = WeightedGraph(Graph(n))

function add_edge!(wg::WeightedGraph, src::Int, dst::Int, weight::Number = 1)
  add_edge!(wg.graph, src, dst)
  push!(wg.weights, (src, dst) => weight)
end

function weights(wg::WeightedGraph)
  return wg.weights
end

nv(wg::WeightedGraph) = nv(wg.graph)

ne(wg::WeightedGraph) = ne(wg.graph)

edges(wg::WeightedGraph) = edges(wg.graph)

induced_subgraph(wg::WeightedGraph, elist) = induced_subgraph(wg.graph, elist)

connected_components(wg::WeightedGraph) = connected_components(wg.graph)
