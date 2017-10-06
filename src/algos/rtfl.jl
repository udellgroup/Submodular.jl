#############################################################################
# rtfl.jl
# Implements the Regularized Follow The Leader algorithm.
# For reference look at:
#############################################################################

export decomposition

function decomposition(w::AbstractExpr, p::BasePoly)
  if !in(w, p)
    error("Cannot do decomposition when the vector is not in the polyhedron.")
  else
    return minimum_norm_point(p, w)[2]
  end
end
