using CombiOpt
using FactCheck

TOL = 1e-3

facts("Prox") do

  context("Associated Polyhedra") do
    # the Lovasz extension of a submodular function defined on a set of two elements
    function twodim(A::AbstractVector)
      if A == []
        return 0
      end
      if A == [1]
        return 3
      end
      if A == [3]
        return 2
      end
      if Set(A) == Set([1, 3])
        return 4
      end
      if in(2, A)
        return 0
      end
    end
    w1 = [3, 1]
    w2 = [2, 2]
    w3 = [0, 0]
    w4 = [3, 2]
    w5 = [0, 5]
    w6 = [-2, -2]
    w7 = [50, -3]
    w8 = [-5, 1]

    p1 = SubPoly(twodim, [1, 3])
    @fact prox(p1, w1) --> roughly([3, 1], TOL)
    @fact prox(p1, w2) --> roughly([2, 2], TOL)
    @fact prox(p1, w3) --> roughly([0, 0], TOL)
    @fact prox(p1, w4) --> roughly([2.5, 1.5], TOL)
    @fact prox(p1, w5) --> roughly([0, 2], TOL)
    @fact prox(p1, w6) --> roughly([-2, -2], TOL)
    @fact prox(p1, w7) --> roughly([3, -3], TOL)
    @fact prox(p1, w8) --> roughly([-5, 1], TOL)

    p2 = BasePoly(twodim, [1, 3])
    @fact prox(p2, w1) --> roughly([3, 1], TOL)
    @fact prox(p2, w2) --> roughly([2, 2], TOL)
    @fact prox(p2, w3) --> roughly([2, 2], TOL)
    @fact prox(p2, w4) --> roughly([2.5, 1.5], TOL)
    @fact prox(p2, w5) --> roughly([2, 2], TOL)
    @fact prox(p2, w6) --> roughly([2, 2], TOL)
    @fact prox(p2, w7) --> roughly([3, 1], TOL)
    @fact prox(p2, w8) --> roughly([2, 2], TOL)

    p3 = PosPoly(twodim, [1, 3])
    @fact prox(p3, w1) --> roughly([3, 1], TOL)
    @fact prox(p3, w2) --> roughly([2, 2], TOL)
    @fact prox(p3, w3) --> roughly([0, 0], TOL)
    @fact prox(p3, w4) --> roughly([2.5, 1.5], TOL)
    @fact prox(p3, w5) --> roughly([0, 2], TOL)
    @fact prox(p3, w6) --> roughly([0, 0], TOL)
    @fact prox(p3, w7) --> roughly([3, 0], TOL)
    @fact prox(p3, w8) --> roughly([0, 1], TOL)

    p4 = SymPoly(twodim, [1, 3])
    @fact prox(p4, w1) --> roughly([3, 1], TOL)
    @fact prox(p4, w2) --> roughly([2, 2], TOL)
    @fact prox(p4, w3) --> roughly([0, 0], TOL)
    @fact prox(p4, w4) --> roughly([2.5, 1.5], TOL)
    @fact prox(p4, w5) --> roughly([0, 2], TOL)
    @fact prox(p4, w6) --> roughly([-2, -2], TOL)
    @fact prox(p4, w7) --> roughly([3, -1], TOL)
    @fact prox(p4, w8) --> roughly([-3, 1], TOL)
  end

end
