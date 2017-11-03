using CombiOpt
using FactCheck

Tol = 1e-3

facts("Algorithms") do

  context("card-inc-fix") do
    n = 20
    g = rand(n)                      # the vector used to generate the cardinality-based function
    g = sort(g, rev=true)
    g = (g + 0.1)/1.01
    for i = 2: n
      g[i] = g[i] + g[i-1]
    end
    S = SetVariable(n)
    gg = card(g, S)                  # F(S) = g(|S|)

    y = rand(n)                      # the point to be projected on the base polytope of F(S) = g(|S|)
    euclidean_proj = card_inc_fix(gg, y, "euclidean")

    # sanity checks:
    # Is the sorted order of indices in y and the projection the same?
    # this is a known property of projections under uniform divergences over cardinality-based polytopes
    @fact sortperm(y) == sortperm(euclidean_proj[1:n]) --> true

    # Is the constraint x(E) = F(E) = g(n) satisfied upto an error of Tol^3?
    @fact sum(euclidean_proj) - g[n] --> roughly(0, Tol^3)
  end

end
