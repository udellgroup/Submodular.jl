#############################################################################
# extreme_point.jl
# Compute the of the base polyhedron given an ordering.
#############################################################################

export extremepoint

function extremepoint(F::Function, w::AbstractArray,
                S = collect(1:length(w))::AbstractArray)
  @assert F([])[1] == 0 # f should be 0 at the empty set
  n = length(w)
  @assert length(S) == n
  x = zeros(n)
  for ii = 1:n
    x[w[ii]] = F(V[w[1: ii]])[1] - F(V[w[1: ii-1]])[1]
  end
  return x
end

extremepoint(F::SubmodFunc, w::AbstractArray) = extremepoint(x -> evaluate(F, x), w)
