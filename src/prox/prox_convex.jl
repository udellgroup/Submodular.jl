#############################################################################
# prox_convex.jl
# compute the prox of convex functions
#############################################################################

export prox

function prox(f::AbstractExpr, w::AbstractArray)
  x = get_cv(f)[1]
  problem = Problem(:minimize, f + 0.5 * norm(x - w)^2)
  solve!(problem, SCSSolver())
  return(x.value)
end
