#############################################################################
# solutions.jl
#############################################################################

export Solution
export solve!

function solve!(p::SCOPEProblem{ConvexProblem}; s::AbstractMathProgSolver = MosekSolver())
  q = Problem(p.head, p.objective, p.constraints)
  solve!(q, s)
  p.status = q.status
  p.optval = q.optval
  # p.solution = q.solution
end

function solve!(p::SCOPEProblem{ConvexLovasz};
                s::AbstractMathProgSolver = MosekSolver(),
                abs_tol::Float64 = 1e-3,
                max_iters::Int = 200,
                max_rep::Int = 100,
                lev_tol_def::Float64 = 1e-3,
                lev_tol_up::Float64 = 1e-1,
                lev_tol_low::Float64 = 1e-7,
                shrink_point::Int = 20,
                λ::Float64 = 1.5,
                μ::Float64 = 0.5)
  cutting_plane(p.model.convex_part, p.model.lovasz, abs_tol, max_iters, max_rep, lev_tol_def, lev_tol_up, lev_tol_low, shrink_point, λ, μ)
end

function solve!(p::SCOPEProblem{ConvexLovaszAbs};
                s::AbstractMathProgSolver = MosekSolver(),
                abs_tol::Float64 = 1e-3,
                max_iters::Int = 200,
                max_rep::Int = 100,
                lev_tol_def::Float64 = 1e-3,
                lev_tol_up::Float64 = 1e-1,
                lev_tol_low::Float64 = 1e-7,
                shrink_point::Int = 20,
                λ::Float64 = 1.5,
                μ::Float64 = 0.5)
  cutting_plane(p.model.convex_part, p.model.lovaszabs, abs_tol, max_iters, max_rep, lev_tol_def, lev_tol_up, lev_tol_low, shrink_point, λ, μ)
end

function solve!(p::SCOPEProblem{LPoverAssocPoly})
  solveLPOA!(p.objective, p.constraints[1].rhs)
  p.optval = evaluate(p.objective)
end

function solve!(p::SCOPEProblem{AssocPolyConstrained};
                            maxiters::Int64 = 100,
  													epsilon::Float64 = 1e-3,
                            verbose::Bool = true)
  frank_wolfe_away(p, maxiters, epsilon, verbose)
end

# function solve!(p::SCOPEProblem;
#                 s::AbstractMathProgSolver = MosekSolver(),
#                 abs_tol::Float64 = 1e-3,
#                 max_iters::Int = 200,
#                 max_rep::Int = 100,
#                 lev_tol_def::Float64 = 1e-3,
#                 lev_tol_up::Float64 = 1e-1,
#                 lev_tol_low::Float64 = 1e-7,
#                 shrink_point::Int = 20,
#                 λ::Float64 = 1.5,
#                 μ::Float64 = 0.5)
#
#   if p.model == ConvexProblem()
#     q = Problem(p.head, p.objective, p.constraints)
#     solve!(q, s)
#     p.status = q.status
#     p.optval = q.optval
#     p.solution = q.solution
#   end
#   if p.model == LPoverAssocPoly()
#     solveLPOA!(p.objective, p.constraints[1].rhs)
#     p.optval = evaluate(p.objective)
#   elseif p.model == AssocPolyConstrained()
#     chambolle_pock(p.objective, p.constraints[1].rhs)
#     p.optval = evaluate(p.objective)
#   elseif isa(p.model, ConvexLovasz)
#     cutting_plane(p.model.convex_part, p.model.lovasz, abs_tol, max_iters, max_rep, lev_tol_def, lev_tol_up, lev_tol_low, shrink_point, λ, μ)
#   elseif isa(p.model, ConvexLovaszAbs)
#     cutting_plane(p.model.convex_part, p.model.lovaszabs, abs_tol, max_iters, max_rep, lev_tol_def, lev_tol_up, lev_tol_low, shrink_point, λ, μ)
#   end
# end
