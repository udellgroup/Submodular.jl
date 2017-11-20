#############################################################################
# solve_dual.jl
# Calculate the dual problem of the given convex programming problems in their
# conic forms.
# Conic problems having the form
# minimize c'*x
# st       b - Ax ∈ K, where K is a cone,
# their dual problems should be
# maximize -b'*y
# st       c + Aᵀy ∈ K^*, where K^* is the dual cone of K.
#############################################################################

import MathProgBase
using Mosek
using ECOS

export solvedualMosek!

function solvedualMosek!(p::Problem)

  c, A, b, dual_var_cones, var_to_ranges, vartypes, conic_constraints = conic_problem(p)
  for i = 1:length(dual_var_cones)
    if dual_var_cones[i][1] == :Free
      dual_var_cones[i] = (:Zero, dual_var_cones[i][2])
    elseif dual_var_cones[i][1] == :Zero
      dual_var_cones[i] = (:Free, dual_var_cones[i][2])
    elseif dual_var_cones[i][1] == :ExpPrimal
      dual_var_cones[i] = (:ExpDual, dual_var_cones[i][2])
    elseif dual_var_cones[i][1] == :ExpDual
      dual_var_cones[i][1] = (:ExpPrimal, dual_var_cones[i][2])
    end
  end
  dual_constr_cones = fill((:Zero, 1:size(A, 2)),1)
  m = ConicModel(MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-6))
  loadproblem!(m, vec(full(b)), -A', vec(full(c)), dual_constr_cones, dual_var_cones)

  optimize!(m)

  solution = try
    MathProgBase.getsolution(m)
  catch
    fill(NaN, numvar(m))
  end

  # return cones, vec(full(c)), A, vec(full(b))
  return dot(vec(full(b)), -solution)
end

# function solvedualECOS!(p::Problem)
#
#   c, A, b, dual_var_cones, var_to_ranges, vartypes, conic_constraints = conic_problem(p)
#   for i = 1:length(dual_var_cones)
#     if dual_var_cones[i][1] == :Free
#       dual_var_cones[i] = (:Zero, dual_var_cones[i][2])
#     elseif dual_var_cones[i][1] == :Zero
#       dual_var_cones[i] = (:Free, dual_var_cones[i][2])
#     elseif dual_var_cones[i][1] == :ExpPrimal
#       dual_var_cones[i] = (:ExpDual, dual_var_cones[i][2])
#     elseif dual_var_cones[i][1] == :ExpDual
#       dual_var_cones[i][1] = (:ExpPrimal, dual_var_cones[i][2])
#     end
#   end
#   dual_constr_cones = fill((:Zero, 1:size(A, 2)),1)
#   m = ConicModel(ECOSSolver(verbose=false, abstol = 1e-4))
#   loadproblem!(m, vec(full(b)), -A', vec(full(c)), dual_constr_cones, dual_var_cones)
#
#   optimize!(m)
#
#   solution = try
#     MathProgBase.getsolution(m)
#   catch
#     fill(NaN, numvar(m))
#   end
#
#   # return cones, vec(full(c)), A, vec(full(b))
#   return dot(vec(full(b)), -solution)
# end
