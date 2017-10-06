#############################################################################
# extreme_point.jl
# Compute the of the base polyhedron given an ordering.
#############################################################################

export extremepoint

function extremepoint(f::Function, w::AbstractArray,
                S = collect(1:length(w))::AbstractArray)
  @assert f([])[1] == 0 # f should be 0 at the empty set
  n = length(w)
  @assert length(S) == n
  x = zeros(n)
  for ii = 1:n
    x[w[ii]] = f(V[w[1: ii]])[1] - f(V[w[1: ii-1]])[1]
  end
  return x
end

extremepoint(f::AbstractExpr, w::AbstractArray) = extremepoint(x -> evaluate(f, x), w)

extremepoint(f::CombiFunc, w::AbstractArray) = extremepoint(x -> evaluate(f, x), w)
