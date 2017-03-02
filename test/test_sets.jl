using CombiOpt
using FactCheck

TOL = 1e-3

facts("Sets") do

  context("Permutation") do
    p = perm(4)
    c = [.4,.2,.6,.3]
    @fact fenchel(p, c) --> roughly([2, 4, 1, 3], TOL)
  end

end
