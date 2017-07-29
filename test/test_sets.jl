using CombiOpt
using FactCheck

TOL = 1e-3

facts("Sets") do

  context("Permutation") do
    p = Perm(4)
    d1 = randperm(4)
    @fact in(d1, p) --> true
    d2 = collect(2:5)
    @fact in(d2, p) --> false
    c = [.4,.2,.6,.3]
    @fact fenchel(p, c) --> roughly([2, 4, 1, 3], TOL)
  end

  context("Polyhedra") do
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
    w1 = [3, 2]
    w2 = [2, 2]
    w3 = [-2, -1]

    p1 = SubPoly(twodim, [3, 1])
    @fact in(w1, p1) --> false
    @fact in(w2, p1) --> true
    @fact in(w3, p1) --> true
    @fact fenchel(p1, w1) --> -Inf
    @fact fenchel(p1, w2) --> -Inf
    @fact fenchel(p1, w3) --> roughly(-7.0, TOL)

    p2 = BasePoly(twodim, [3, 1])
    @fact in(w1, p2) --> false
    @fact in(w2, p2) --> true
    @fact in(w3, p2) --> false
    @fact fenchel(p2, w1) --> roughly(10.0, TOL)
    @fact fenchel(p2, w2) --> roughly(8.0, TOL)
    @fact fenchel(p2, w3) --> roughly(-7.0, TOL)

    p3 = PosPoly(twodim, [3, 1])
    @fact in(w1, p3) --> false
    @fact in(w2, p3) --> true
    @fact in(w3, p3) --> false
    @fact fenchel(p3, w1) --> roughly(0.0, TOL)
    @fact fenchel(p3, w2) --> roughly(0.0, TOL)
    @fact fenchel(p3, w3) --> roughly(-7.0, TOL)

    p4 = SymPoly(twodim, [3, 1])
    @fact in(w1, p4) --> false
    @fact in(w2, p4) --> true
    @fact in(w3, p4) --> true
    @fact fenchel(p4, w1) --> roughly(-11.0, TOL)
    @fact fenchel(p4, w2) --> roughly(-8.0, TOL)
    @fact fenchel(p4, w3) --> roughly(-7.0, TOL)
  end

end
