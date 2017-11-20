#############################################################################
# weights_of_graphs.jl
# handles the weights of AbstractGraph, complementary of LightGraphs.
#############################################################################

export WeightedGraph
export add_edge!, vertice_weights, edge_weights
export nv, ne, edges, induced_subgraph, connected_components

type WeightedGraph
  vertice_weights::Dict{Int64, Number}
  edge_weights::Dict{Tuple{Int64,Int64},Number}
  graph::AbstractGraph

  function WeightedGraph(graph::AbstractGraph, vertice_weights::Dict{Int64, Number} = Dict{Int64, Number}())
    for i = 1:nv(graph)
      if in(i, keys(vertice_weights)) == false
        push!(vertice_weights, i => 1)
      end
    end
    return new(vertice_weights, Dict{Tuple{Int64,Int64}, Number}(), graph)
  end
end

WeightedGraph(n::Int) = WeightedGraph(Graph(n))

WeightedGraph(n::Int, vertice_weights::Dict{Int64, Number}) = WeightedGraph(Graph(n), vertice_weights)

function add_edge!(wg::WeightedGraph, src::Int, dst::Int, weight::Number = 1)
  add_edge!(wg.graph, src, dst)
  push!(wg.edge_weights, (src, dst) => weight)
end

function vertice_weights(wg::WeightedGraph)
  return wg.vertice_weights
end

function edge_weights(wg::WeightedGraph)
  return wg.edge_weights
end

nv(wg::WeightedGraph) = nv(wg.graph)

ne(wg::WeightedGraph) = ne(wg.graph)

edges(wg::WeightedGraph) = edges(wg.graph)

induced_subgraph(wg::WeightedGraph, elist) = induced_subgraph(wg.graph, elist)

connected_components(wg::WeightedGraph) = connected_components(wg.graph)
