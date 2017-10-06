#############################################################################
# solvers.jl
#############################################################################

export SCOPESolver

abstract type SCOPESolver end

type ConvexProblem <: SCOPESolver end

# Solve the problem as a convex programming problem if its structure is not identified
function get_solver(objective::AbstractExpr, constraints::Array{Constraint})
  return ConvexProblem()
end
