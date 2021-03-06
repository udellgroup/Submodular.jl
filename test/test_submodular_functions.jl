using Submodular, LightGraphs
using FactCheck

Tol = 1e-3

facts("Submodular Functions") do

  context("card-based atom") do
    S = SetVariable(3)
    F = card(S)
    @fact evaluate(F, [1, 3]) --> roughly(2, Tol)

    p(z) = -0.5*z^2 + 3*z + 0.5 * z
    perm_func = compose(p, F)
    @fact modularity(perm_func) --> SubModularity()
    @fact evaluate(perm_func, [1, 2, 3])[1] --> roughly(6, Tol)

    F₁ = card([3.0, 5, 6], S)
    @fact modularity(F₁) --> SubModularity()
    @fact evaluate(F₁, [1, 3]) --> roughly(5, Tol)
  end

  context("customized submodular function atom") do
    function func(S::AbstractArray)
      if S == []
        return 0
      elseif S == [1]
        return 3
      elseif S == [2]
        return 2
      elseif S == [1, 2]
        return 4
      end
    end
    S = SetVariable(2)
    F = submod(func, S)
    @fact modularity(F) --> SubModularity()
    @fact evaluate(F, [1]) --> roughly(3, Tol)
    @fact evaluate(F, [1, 2]) --> roughly(4, Tol)
  end

  context("cut atom and WeightedGraph") do
    S = SetVariable(3)
    G = WeightedGraph(3)
    add_edge!(G, 1, 3)
    F = cut(G, S)
    @fact modularity(F) --> SubModularity()
    @fact evaluate(F, [1, 2]) --> roughly(1, Tol)
  end

  context("log determinant") do
    S = SetVariable(3)
    M = [1 0 0; 0 2 0; 0 0 3]
    F = logdet(M, S)
    @fact modularity(F) --> SubModularity()
    @fact evaluate(F, [1, 3]) --> roughly(log(3), Tol)
  end

  context("modular atom") do
    S = SetVariable(4)
    F = modular([1, 2, 3, 4], S)
    @fact modularity(F) --> Modularity()
    @fact evaluate(F, [1, 3]) --> roughly(4, Tol)
  end

  context("rank of graph matroid atom and WeightedGraph") do
    S = SetVariable(2)
    G = WeightedGraph(3)
    add_edge!(G, 1, 3)
    add_edge!(G, 2, 3)
    F = rank_of_graph_matroid(G, S)
    @fact modularity(F) --> SubModularity()
    @fact evaluate(F, [1, 2]) --> roughly(2, Tol)
  end

end
