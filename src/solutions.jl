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
                solver::String = "cutting_plane",
                s::AbstractMathProgSolver = MosekSolver(),
                abs_tol::Float64 = 1e-3,
                max_iters::Int = 100,
                max_rep::Int = 100,
                lev_tol_def::Float64 = 1e-3,
                lev_tol_up::Float64 = 1e-1,
                lev_tol_low::Float64 = 1e-7,
                shrink_point::Int = 20,
                λ::Float64 = 1.5,
                μ::Float64 = 0.7,              # linesearch shrink rate or
                x₀::AbstractArray = zeros(length(get_v(p.objective)[1])),   # the starting point for L_BFGS
                m::Int = length(get_v(p.objective)[1]),  # the number of kept iterations
                c₁::Float64 = 0.5,             # linesearch parameter
                c₂::Float64 = 0.1,             # linesearch parameter
                epsilon::Float64 = 1e-3)       # the stopping criteria of the algorithm

  if solver == "cutting_plane"
    cutting_plane(p.model.convex_part, p.model.lovasz, abs_tol = abs_tol, max_iters = max_iters, max_rep = max_rep, lev_tol_def = lev_tol_def, lev_tol_up = lev_tol_up, lev_tol_low = lev_tol_low, shrink_point = shrink_point, λ = λ, μ = μ)
  elseif solver == "L_BFGS"
    L_BFGS(p.model.convex_part, p.model.lovasz, m = m, c₁ = c₁, c₂ = c₂, μ = μ, epsilon = epsilon, max_iters = max_iters)
  else
    error("This method is not supported.")
  end
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
                solver::String = "frank_wolfe_away",
                max_iters::Int64 = 100,
  							epsilon::Float64 = 1e-3,
                verbose::Bool = true,
                s::AbstractMathProgSolver = MosekSolver(),
                abs_tol::Float64 = 1e-3,
                max_rep::Int = 100,
                lev_tol_def::Float64 = 1e-3,
                lev_tol_up::Float64 = 1e-1,
                lev_tol_low::Float64 = 1e-7,
                shrink_point::Int = 20,
                λ::Float64 = 1.5,
                μ::Float64 = 0.7)

  if solver == "frank_wolfe_away"
    frank_wolfe_away(p, max_iters = max_iters, epsilon = epsilon, verbose = verbose)
  elseif solver == "cutting_plane"
    cutting_plane_dual
  else
    error("This method is not supported.")
  end
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
