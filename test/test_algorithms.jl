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

  context("frank-wolfe with away steps") do
    n = 4
    x = Variable(n)
    g = norm(x - [3, 2, 5, 1])
    S = SetVariable(n)
    F = card(S)
    p(z) = -0.5*z^2 + n*z + 0.5 * z
    perm_func = compose(p, F)
    P = BasePoly(perm_func)
    prob = SCOPEminimize(g, x in P)
    @fact frank_wolfe_away(prob, verbose = false) --> roughly([3, 2, 4, 1], Tol)
  end

end
