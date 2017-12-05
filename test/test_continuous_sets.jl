using Submodular
using FactCheck

Tol = 1e-3

facts("Continuous Sets") do

  context("permutation atom") do
    p = Perm(4)
    d1 = randperm(4)
    @fact in(d1, p) --> true
    d2 = collect(2:5)
    @fact in(d2, p) --> false
    w = [.4,.2,.6,.3]
    @fact fenchel(p, w) --> roughly([2, 4, 1, 3], Tol)
  end

  context("associated polyhedra and the Fujishige-minimum-norm-point algorithm") do
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
    @fact in([1, 2, 5, 3], P₂) --> false
    @fact in([1, 2, 4, 3], P₂) --> true
    @fact in([1, 2, 3, 3], P₂) --> false
    @fact in([-1, -2, 0, 0], P₂) --> false
    @fact in([-2, -2, -5, -3], P₂) --> false

    P₃ = PosPoly(perm_func)
    @fact in([1, 2, 5, 3], P₃) --> false
    @fact in([1, 2, 4, 3], P₃) --> true
    @fact in([1, 2, 3, 3], P₃) --> true
    @fact in([-1, -2, 0, 0], P₃) --> false
    @fact in([-2, -2, -5, -3], P₃) --> false

    P₄ = SymPoly(perm_func)
    @fact in([1, 2, 5, 3], P₄) --> false
    @fact in([1, 2, 4, 3], P₄) --> true
    @fact in([1, 2, 3, 3], P₄) --> true
    @fact in([-1, -2, 0, 0], P₄) --> true
    @fact in([-2, -2, -5, -3], P₄) --> false
  end

end
