import Base.in

export SetConstraint
export in

### Set constraint
mutable struct SetConstraint{T <: ContiSet} <: Constraint
  head::Symbol
  id_hash::UInt64
  lhs::AbstractExpr
  rhs::T
  size::Tuple{Int, Int}
  dual::ValOrNothing
end

function SetConstraint(lhs::Variable, rhs::ContiSet)
  id_hash = hash((lhs, rhs, :(==)))
  return SetConstraint(:(set), id_hash, lhs, rhs, (1, 1), nothing)
end

function get_cv(constraint::SetConstraint)
  return get_cv(constraint.lhs)
end

function get_sv(constraint::SetConstraint)
  return get_sv(constraint.lhs)
end

in(lhs::Variable, rhs::ContiSet) = SetConstraint(lhs, rhs)

# the fenchel conjugate of a set indicator function
fenchel(setconstraint::SetConstraint, w::AbstractArray) = fenchel(setconstraint.rhs, w)

# the prox of a set indicator function
prox(setconstraint::SetConstraint, w::AbstractArray) = prox(setconstraint.rhs, w)
