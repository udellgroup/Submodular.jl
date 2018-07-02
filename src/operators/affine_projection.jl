#############################################################################
# aff_projection.jl
# Projects w onto the affine hull of point in S
# Note that points in S should be linearly independent
#############################################################################

export aff_proj

function aff_proj(w::AbstractArray, S::AbstractMatrix)
  T = S[:, 2: end] - S[:, 1] * ones(1, size(S, 2) - 1)
  x1 = S[:, 1]
  y = T*((T'*T)\T')*(w - x1) + x1
  return y
end
