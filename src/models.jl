#############################################################################
# models.jl
# Handles the models of problems
#############################################################################

export CombiModel, LPoverAssocPoly, ConvexLovasz, AssocPolyConstrained

function get_model(objective::AbstractExpr, constraints::AbstractArray=Constraint[])
  objective_cv = get_cv(objective)
  objective_sv = get_sv(objective)
  constraints_cv = []
  constraints_sv = []
  for i = 1:length(constraints)
    constraints_cv = vcat(constraints_cv, get_cv(constraints[i].lhs), get_cv(constraints[i].rhs))
    constraints_sv = vcat(constraints_sv, get_sv(constraints[i].lhs), get_sv(constraints[i].rhs))
  end
  if objective.head == :+
    lonum = 0
    loind = 1
    lanum = 0
    laind = 1
    childnum = length(objective.children)
    for i = 1:childnum
      if objective.children[i].head == :lovasz
        lonum += 1
        loind = i
      elseif objective.children[i].head == :lovaszabs
        lanum += 1
        laind = i
      end
    end
    if (lonum + lanum) > 1
      error("Cannot solve objective functions with multiple Lovasz extensions!")
    elseif lonum == 1
      lovasz = objective.children[loind]
      if loind == 1
        objective1 = objective.children[2]
        if childnum > 2
          for i = 3 : childnum
            objective1 += objective.children[i]
          end
        end
      else
        objective1 = objective.children[1]
        for i = 2:(loind - 1)
          objective1 += objective.children[i]
        end
        for i = (loind + 1) : childnum
          objective1 += objective.children[i]
        end
      end
      if vexity(objective1) == ConvexVexity() || vexity(objective1) == AffineVexity()
        q = Problem(:minimize, objective1, constraints)
        return ConvexLovasz(q, lovasz)
      else
        warn("The objective function is the sum of the Lovasz extension of a submodular function and a non-DCP-convex funcion; solution may be inaccurate.")
        q = Problem(:minimize, objective1, constraints)
        return ConvexLovasz(q, lovasz)
      end
    elseif lanum == 1
      lovasz = objective.children[laind]
      if laind == 1
        objective1 = objective.children[2]
        if childnum > 2
          for i = 3 : childnum
            objective1 += objective.children[i]
          end
        end
      else
        objective1 = objective.children[1]
        for i = 2:(laind - 1)
          objective1 += objective.children[i]
        end
        for i = (laind + 1) : childnum
          objective1 += objective.children[i]
        end
      end
      if vexity(objective1) == ConvexVexity() || vexity(objective1) == AffineVexity()
        q = Problem(:minimize, objective1, constraints)
        return ConvexLovaszAbs(q, lovasz)
      else
        warn("The objective function is the sum of the Lovasz extension of a submodular function and a non-DCP-convex funcion; solution may be inaccurate.")
        q = Problem(:minimize, objective1, constraints)
        return ConvexLovaszAbs(q, lovasz)
      end
    else
      return ConvexProblem()
    end
  elseif length(constraints) == 1
    if typeof(constraints[1]) <: SetConstraint
      if typeof(constraints[1].rhs) <: AssocPoly && length(objective_sv) == 0 && length(objective_cv) == 1 && objective_cv[1].id_hash == constraints_cv[1].id_hash
        # if typeof(vexity(objective)) <: AffineVexity
        if vexity(objective) == AffineVexity()
          return LPoverAssocPoly()
        else
          return AssocPolyConstrained()
        end
      end
    else
      return ConvexProblem()
    end
  else
    return ConvexProblem()
  end
end

### convex problem
type ConvexProblem <: SCOPEModel end
### convex optimization over a polyhedron associated with a submodular function
type AssocPolyConstrained <: SCOPEModel end
### linear programming over associated polyhedra
type LPoverAssocPoly <: SCOPEModel end
