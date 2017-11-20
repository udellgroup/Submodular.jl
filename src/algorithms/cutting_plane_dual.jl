#############################################################################
# cutting_plane_dual.jl
# convert a problem of the DualModel to the PrimalModel and solve it with the
# cutting plane method.
#############################################################################

export cutting_plane_dual

function cutting_plane_dual(p::SCOPEProblem{AssocPolyConstrained},
                            f::LovaszExtAtom;
                            s::AbstractMathProgSolver = MosekSolver(),
                            abs_tol::Float64 = 1e-3,
                            max_iters::Int = 100,
                            max_rep::Int = 100,
                            lev_tol_def::Float64 = 1e-3,
                            lev_tol_up::Float64 = 1e-1,
                            lev_tol_low::Float64 = 1e-7,
                            shrink_point::Int = 20,
                            λ::Float64 = 1.5,
                            μ::Float64 = 0.7)
end

function dual_to_primal(p::Problem)
  c, A, b, dual_var_cones, var_to_ranges, vartypes, conic_constraints = conic_problem(p)
  # convert the problem into its dual form
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
  c, A, b = vec(full(b)), -A', vec(full(c))
  push!(c, 1)
