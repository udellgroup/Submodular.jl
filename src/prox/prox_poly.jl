#############################################################################
# prox_poly.jl
# compute the prox of indicator functions of associated polyhedra
#############################################################################

export prox

prox(p::AssocPoly, w::AbstractArray) = minimum_norm_point(p, w)

prox(C::SetConstraint{AssocPoly}, w::AbstractArray) = minimum_norm_point(p, w)
