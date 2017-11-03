using CombiOpt
using Convex
using FactCheck

Tol = 1e-3

facts("Operators") do

  context("affine projection") do
    S = [3 2; 1 2]
    w1 = [2, 2]
    w2 = [3, 2]
    w3 = [5, 0]
    @fact affproj(w1, S) --> roughly([2, 2], Tol)
    @fact affproj(w2, S) --> roughly([2.5, 1.5], Tol)
    @fact affproj(w3, S) --> roughly([4.5, -0.5], Tol)
  end

  context("extreme point") do
    S = SetVariable(4)
    F = card([40, 70, 90, 100], S)
    @fact extreme_point(F, [3, 2, 4, 1]) --> roughly([10, 30, 40, 20], Tol)
  end

  context("gradient atom") do
    x = Variable(2)
    f = norm(x)
    g = grad(f)
    a = rand(2)
    @fact evaluate(g, a) * norm(a) --> roughly(a, Tol)
  end

  context("greedy algorithm") do
    S = SetVariable(4)
    F = card([40, 70, 90, 100], S)
    @fact greedy(F, [35, 25, 45, 15]) --> roughly([30, 20, 40, 10], Tol)
  end

  context("Lovasz extension") do
    # the Lovasz extension of a submodular function defined on a set of two elements
    S = SetVariable(4)
    F = modular([1, 2, 3, 4], S)
    x = Variable(4)
    f = lovasz(F, x)
    @fact vexity(f) --> AffineVexity()
    @fact evaluate(f, [1.5, 3.5, 2.5, 4.5]) --> roughly(34, Tol)
  end

  context("Lovasz extension of absolute value") do
    # the Lovasz extension of a submodular function defined on a set of two elements
    S = SetVariable(4)
    F = modular([1, 2, 3, 4], S)
    x = Variable(4)
    f = lovasz(F, abs(x))
    @fact vexity(f) --> ConvexVexity()
    @fact evaluate(f, [-1.5, 3.5, -2.5, 4.5]) --> roughly(34, Tol)
  end

end
