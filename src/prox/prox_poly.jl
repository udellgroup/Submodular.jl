#############################################################################
# prox_poly.jl
# compute the prox of indicator functions of associated polyhedra
#############################################################################

export prox

function prox(p::AssocPoly, w::AbstractArray, Tol::Float64 = 1e-3)
  if isa(p, BasePoly)
    if isa(p.F, CardBasedAtom)
      return  card_inc_fix(p.F, w)
    end
  else
    minimum_norm_point(p, w, Tol)
  end
end

prox(C::SetConstraint{AssocPoly}, w::AbstractArray, Tol::Float64 = 1e-3) = prox(p, w, Tol)
