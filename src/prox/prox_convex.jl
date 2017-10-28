#############################################################################
# prox_convex.jl
# compute the prox of convex functions
#############################################################################

using Mosek

export prox

function prox(f::AbstractExpr, w::AbstractArray, abs_tol::Float64 = 1e-3)
  x = get_cv(f)[1]
  problem = Problem(:minimize, f + 0.5 * norm(x - w)^2)
  solve!(problem, MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = abs_tol))
  return(x.value)
end
