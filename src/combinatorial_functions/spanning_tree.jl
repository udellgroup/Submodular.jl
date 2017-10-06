#############################################################################
# spanning_tree
# Handles the spanning tree function: f: 2ⱽ → ℝ, f(S) = n - connected_components(S),
# where underlying graphs (V, E) are given, and |V| = n.
#############################################################################

export spanning_tree
export sign, monotonicity, modularity, evaluate

type SpanTreeAtom <: CombiFunc   # TODO: extend to sets other than plain set varables
  head::Symbol
  id_hash::UInt64
  children::Tuple{CombiSet}
  size::Tuple{Int, Int}
  param::AbstractExpr
  setvariables::Array{CombiSet}
  graph::AbstractGraph

  function SpanTreeAtom(graph::AbstractGraph, S::CombiSet)
    if collect(g.vertices) != S.baseset
      error("Cannot define a spanning tree function when vertices of the graph and the base set are different.")
    else
      children = (S,)
      setvariables = get_sv(S)
      return new(:spaningtree, hash(children), children, (1, 1), x, setvariables)
    end
  end
end

spanning_tree(graph::AbstractGraph, S::CombiSet) = SpanTreeAtom(g, S)

function sign(x::SpanTreeAtom)
  return Positive()
end

function monotonicity(x::SpanTreeAtom)
  return NoMonotonicity()
end

function modularity(x::SpanTreeAtom)
  return Nonincreasing()
end

function evaluate(x::SpanTreeAtom) # TODO only allow scalar input
  return x.children[1].cardinality - connected_components
end
