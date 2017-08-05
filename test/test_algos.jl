using CombiOpt
using FactCheck

TOL = 1e-3

facts("Algorithms") do

  context("Minimum Norm") do
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
    V = [1, 3]
    w1 = [3, 1]
    w2 = [2, 2]
    w3 = [0, 0]
    w4 = [3, 2]
    w5 = [0, 5]
    w6 = [-2, -2]
    w7 = [50, -3]
    w8 = [-5, 1]

    @fact minimum_norm(twodim, V, w1) --> roughly([3, 1], TOL)
    @fact minimum_norm(twodim, V, w2) --> roughly([2, 2], TOL)
    @fact minimum_norm(twodim, V, w3) --> roughly([2, 2], TOL)
    @fact minimum_norm(twodim, V, w4) --> roughly([2.5, 1.5], TOL)
    @fact minimum_norm(twodim, V, w5) --> roughly([2, 2], TOL)
    @fact minimum_norm(twodim, V, w6) --> roughly([2, 2], TOL)
    @fact minimum_norm(twodim, V, w7) --> roughly([3, 1], TOL)
    @fact minimum_norm(twodim, V, w8) --> roughly([2, 2], TOL)
  end

end
