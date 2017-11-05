#############################################################################
# prox_poly.jl
# compute the prox of indicator functions of associated polyhedra
#############################################################################

export prox

prox(p::AssocPoly, w::AbstractArray, Tol::Float64 = 1e-3) = minimum_norm_point(p, w, Tol)

prox(C::SetConstraint{AssocPoly}, w::AbstractArray, Tol::Float64 = 1e-3) = minimum_norm_point(p, w, Tol)
