#############################################################################
# extreme_point.jl
# Compute the of the base polyhedron for an given ordering, i.e., the elements
# of the result is calculated in the order of the given ordering.
#############################################################################

export extreme_point

function extreme_point(F::Function, w::AbstractArray,
                       V = collect(1:length(w))::AbstractArray)
  @assert F([])[1] == 0 # f should be 0 at the empty set
  w = Int.(w)
  n = length(w)
  @assert length(V) == n
  x = zeros(n)
  for ii = 1:n
    x[V[w[ii]]] = F(V[w[1: ii]])[1] - F(V[w[1: ii-1]])[1]
  end
  return x
end

extreme_point(F::SubmodFunc, w::AbstractArray) = extreme_point(x -> evaluate(F, x), w)
