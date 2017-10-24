using CombiOpt
using FactCheck

TOL = 1e-3

facts("Continuous Sets") do

  context("permutation atom") do
    p = Perm(4)
    d1 = randperm(4)
    @fact in(d1, p) --> true
    d2 = collect(2:5)
    @fact in(d2, p) --> false
    w = [.4,.2,.6,.3]
    @fact fenchel(p, w) --> roughly([2, 4, 1, 3], TOL)
  end

  context("associated polyhedra") do
    S = SetVariable(4)
    p(z) = -0.5*z^2 + 4*z + 0.5 * z
    perm_func = compose(p, card(S))

    P₁ = SubmodPoly(perm_func)
    @fact in([1, 2, 5, 3], P₁) --> false
    @fact in([1, 2, 4, 3], P₁) --> true
    @fact in([1, 2, 3, 3], P₁) --> true
    @fact in([-1, -2, 0, 0], P₁) --> true
    @fact in([-2, -2, -5, -3], P₁) --> true

    P₂ = BasePoly(perm_func)
    @fact prox(P₂, [1, 2, 5, 3]) --> roughly([1, 2, 4, 3], Tol)
    @fact prox(P₂, [1, 2, 4, 3]) --> roughly([1, 2, 4, 3], Tol)
    @fact prox(P₂, [1, 2, 3, 3]) --> roughly([1.25, 2.25, 3.25, 3.25], Tol)
    @fact prox(P₂, [-1, -2, 0, 0]) --> roughly([2.25, 1.25, 3.25, 3.25], Tol)
    @fact prox(P₂, [-2, -2, -5, -3]) --> roughly([3.33333, 3.33333, 1, 2.33333], Tol)
    #
    # P₃ = PosPoly(perm_func)
    # @fact prox(P₂, [1, 2, 5, 3]) --> roughly([1, 2, 4, 3], Tol)
    # @fact prox(P₂, [1, 2, 4, 3]) --> roughly([1, 2, 4, 3], Tol)
    # @fact prox(P₂, [1, 2, 3, 3]) --> roughly([1.25, 2.25, 3.25, 3.25], Tol)
    # @fact prox(P₂, [-1, -2, 0, 0]) --> roughly([0, 0, 0, 0], Tol)
    # @fact prox(P₂, [-2, -2, -5, -3]) --> roughly([0, 0, 0, 0], Tol)
    #
    # P₄ = SymPoly(perm_func)
    # @fact in([1, 2, 5, 3], P₄) --> false
    # @fact in([1, 2, 3, 3], P₄) --> true
    # @fact in([1, 2, 3, 3], P₄) --> true
    # @fact in([-1, -2, 0, 0], P₄) --> true
    # @fact in([2, 2, 5, 3], P₄) --> false
  end

end
