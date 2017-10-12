#############################################################################
# solutions.jl
#############################################################################

export Solution
export solve!, solve1!

function solve!(p::SCOPEPrimal)
  if p.model == ConvexProblem()
    q = Problem(p.head, p.objective, p.constraints)
    solve!(q)
    p.status = q.status
    p.optval = q.optval
    p.solution = q.solution
  end

  if p.model == LPoverAssocPoly()
    solveLPOA!(p.objective, p.constraints[1].rhs)
    p.optval = evaluate(p.objective)
  elseif p.model == AssocPolyConstrained()
    chambolle_pock(p.objective, p.constraints[1].rhs)
    p.optval = evaluate(p.objective)
  elseif isa(p.model, ConvexLovasz)
    # cutting_plane(p.model)
    cutting_plane(p.model.convex_part, p.model.lovasz)
  elseif isa(p.model, ConvexLovaszAbs)
    cutting_plane(p.model.convex_part, p.model.lovaszabs)
  end
end
