#############################################################################
# lp_over_assoccpoly.jl
# solves linear programming over associated polyhedra.
#############################################################################

export solveLPOA!

function solveLPOA!(objective::AbstractExpr, p::AssocPoly)
  if objective.head == :*
    if isa(objective.children[2], Variable)
      if isa(objective.children[1], Constant)
        objective.children[2].value = fenchel(p, objective.children[1].value[1, :])
      else
        objective.children[2].value = fenchel(p, evaluate(objective.children[1]))
      end
      return objective.children[2].value
    elseif isa(objective.children[1], Variable)
      if isa(objective.children[2], Constant)
        objective.children[1].value = fenchel(p, objective.children[2].value[1, :])
      else
        objective.children[1].value = fenchel(p, evaluate(objective.children[2]))
      end
      return objective.children[1].value
    end
  end
end
