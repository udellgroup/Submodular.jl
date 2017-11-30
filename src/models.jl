#############################################################################
# models.jl
# Handles the models of problems
#
# a PrimalModel has the form:
#   minimize g(x) + f(x)
#   s.t.     x \in s, s in ContiSet
#            x \in C, C is convex,
# where g(x) is a convex function, and f(x) is the Lovasz extension of a submodular function
#
# a DualModel has the form:
#  minimize g(x)
#  s.t.      x \in P, P is a polyhedron associated with a submodular function
##            x \in C, C is convex
#############################################################################

export PrimalModel
export ConvexProblem, ConvexLovasz, ConvexLovaszAbs
export DualModel
export AssocPolyConstrained, LPoverAssocPoly
export get_model

### models of the primal problems
abstract type PrimalModel <: SCOPEModel end

### convex problem
type ConvexProblem <: PrimalModel end
### Convex function + Lovasz extension
type ConvexLovasz <: PrimalModel
  convex_part::Problem
  lovasz::LovaszExtAtom
end
### Convex function + Lovasz extension on absolute values
type ConvexLovaszAbs <: PrimalModel
  convex_part::Problem
  lovaszabs::LovaszExtAbsAtom
end

### models of the dual problems
abstract type DualModel <: SCOPEModel end

### convex optimization over a polyhedron associated with a submodular function
type AssocPolyConstrained <: DualModel
  prob::Problem
  poly::AssocPoly
end
### linear programming over associated polyhedra
type LPoverAssocPoly <: DualModel end

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
    lonum = 0                                         # number of lovasz extensions
    loind = 1                                         # index of lovasz extensions
    lanum = 0                                         # number of lovasz + abs
    laind = 1                                         # index of lovasz + abs
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
      # convex + Lovasz extension
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
      # convex + Lovasz extension on absolute values
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
    if isa(constraints[1], SetConstraint)             # set-constrained optimization problems
      if isa(constraints[1].rhs, AssocPoly) && length(objective_sv) == 0 && length(objective_cv) == 1 && objective_cv[1].id_hash == constraints_cv[1].id_hash
        if vexity(objective) == AffineVexity()
          return LPoverAssocPoly()
        elseif vexity(objective) == ConvexVexity()
          return AssocPolyConstrained(Problem(:minimize, objective), constraints[1].rhs)
        end
      end
    else
      return ConvexProblem()
    end
  else
    return ConvexProblem()                            # fail-safe option, TODO: add other one
  end
end
