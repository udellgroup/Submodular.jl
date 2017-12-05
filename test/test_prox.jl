using Submodular, Convex
using FactCheck

Tol = 1e-3

facts("Prox") do

  context("prox of associated polyhedra") do
    S = SetVariable(4)
    p(z) = -0.5*z^2 + 4*z + 0.5 * z
    perm_func = compose(p, card(S))
    P = SubmodPoly(perm_func)

    x = Variable(4)
    w₁ = [1, 2, 5, 3]
    w₂ = [-2, -2, -5, -3]
    @fact prox(x in P, w₁) --> roughly([1, 2, 4, 3], Tol)
    @fact prox(x in P, w₂) --> roughly([-2, -2, -5, -3], Tol)
  end

  context("prox of convex functions") do
    w = rand(4)
    x = Variable(4)
    f = 0.5 * norm(x)^2
    prox(f, w)
    @fact evaluate(0.5 * w - x) --> roughly(zeros(4, 1), Tol)
  end

end
