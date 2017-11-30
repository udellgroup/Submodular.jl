#############################################################################
# chambolle_pock.jl
# Uses the chambolle_pock algorithm to solve optimization problems consisting
# of two proxable parts f and g.
# For reference look at:
#
# function combi_split(problem, params = CombiParams())
#   @assert ...
#   combi_constraint = problem.combi_constraints[1]
#   z = zeros(size(problem.variable))
#   for t = 1:params.max_iters
#     # compute proximal operator of objective at z
#     reg_objective = problem.objective + 1/2*sumsquares(problem.variable - z)
#     conv_prob = minimize(reg_objective, problem.convex_constraints) # = prox_{objective}(z)
#     solve!(conv_prob)
#     x = evaluate(problem.variable)
#
#     # compute proximal operator of combi_constraint at x
#     z = prox!(combi_constraint, x)
#   end
#   fix!(problem.variable, z) # does this work?
#   problem.solution = z
#   return z
# end
#############################################################################

TOL = 1e-3

export chambolle_pock

function chambolle_pock(f::AbstractExpr, p::AssocPoly,
                        max_iters = 10000 :: Number)
  # Initialization
  x = get_cv(f)[1]
  x0 = zeros(x.size[1])
  u0 = zeros(x.size[1])
  para = Variable(x.size[1])

  # First iteration, in preparation of warm starts
  fix!(para, x0 - u0)
  convex_prob = Problem(:minimize, f + 0.5 * norm(x - para)^2)
  solve!(convex_prob, SCSSolver(verbose = false))
  x1 = (x.value)[:, 1]
  x0 = copy(x1)
  u0 = (u0 + 2 * x1 - x0) - prox(p, (u0 + 2 * x1 - x0), u0)

  # Step 2 iterations
  for iter = 1:max_iters
    fix!(para, x0 - u0)
    solve!(convex_prob, SCSSolver(verbose=false), warmstart = true)
    x1 = (x.value)[:, 1]
    u0 = (u0 + 2 * x1 - x0) - prox(p, (u0 + 2 * x1 - x0), u0)
    if sum(abs.(x1 - x0)) < TOL
      break
    end
    x0 = copy(x1)
  end
  x.value = x0
  return x0
end
