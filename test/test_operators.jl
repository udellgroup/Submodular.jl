using CombiOpt
using FactCheck

TOL = 1e-3

facts("Operators") do

  context("Lovasz Extension") do
    # the Lovasz extension of a submodular function defined on a set of two elements
    function twodim(A::AbstractArray)
      if A == []
        return 0
      end
      if A == [1]
        return 3
      end
      if A == [2]
        return 2
      end
      if Set(A) == Set([1, 2])
        return 4
      end
    end
    w1 = [.5, 1]
    w2 = [1, .5]
    @fact lovasz(twodim, w1) - dot([2, 2], w1) --> roughly(0, TOL)
    @fact lovasz(twodim, w2) - dot([3, 1], w2) --> roughly(0, TOL)

    # the Lovasz extionsion of a modular function is a linear function
    c = rand(4)
    function modular(V::AbstractArray)
      if length(V) == 0
        return 0
      else
        return sum(c[V])
      end
    end
    A = [1, 2, 4]
    x = rand(3)
    @fact lovasz(modular, x, A) - dot(c[A], x)  --> roughly(0, TOL)
  end

  context("Affine Projection") do
    S = [3 2; 1 2]
    w1 = [2, 2]
    w2 = [3, 2]
    w3 = [5, 0]
    @fact affpro(w1, S) --> roughly([2, 2], TOL)
    @fact affpro(w2, S) --> roughly([2.5, 1.5], TOL)
    @fact affpro(w3, S) --> roughly([4.5, -0.5], TOL)
  end

end
