using CombiOpt
using FactCheck

Tol = 1e-3

facts("Combinatorial Sets") do

  context("intersection") do
    A = SetVariable(2)
    B = SetVariable([1, 4, 5])
    C = [1, 10, 11]
    fix!(A, [1, 2])
    fix!(B, [1, 4])
    D = intersect(A, B, C)
    @fact Set(get_elements(D)) == Set([1]) --> true
    @fact Set(D.baseset) == Set([1, 2, 4, 5]) --> true
  end

  context("set difference") do
    A = SetVariable(2)
    B = SetVariable([1, 4, 5])
    fix!(A, [1, 2])
    fix!(B, [1, 4])
    C = setdiff(A, B)
    @fact Set(get_elements(C)) == Set([2]) --> true
    @fact Set(C.baseset) == Set([1, 2, 4, 5]) --> true


    B₁ = setdiff(A, [1, 3])
    @fact Set(get_elements(B₁)) == Set([2]) --> true
    @fact Set(B₁.baseset) == Set([1, 2]) --> true

    B₂ = setdiff([1, 3], A)
    @fact Set(get_elements(B₂)) == Set([3]) --> true
    @fact Set(B₂.baseset) == Set([1, 2]) --> true
  end

  context("union") do
    A = SetVariable(3)
    fix!(A, [1, 2])
    B = union(A, [5])
    @fact Set(get_elements(B)) == Set([1, 2, 5]) --> true
    @fact Set(B.baseset) == Set([1, 2, 3, 5]) --> true

    A₁ = SetVariable(5)
    fix!(A₁, [5])
    B₁ = union(A, A₁)
    @fact Set(get_elements(B₁)) == Set([1, 2, 5]) --> true
    @fact Set(B₁.baseset) == Set([1, 2, 3, 4, 5]) --> true
  end

end
