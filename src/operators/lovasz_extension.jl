# compute the Lovasz extension of f, using the greedy algorithm

export lovaszext

function lovaszext(f::Function, w::AbstractVector,
                   S = collect(1:length(w))::AbstractVector)
  # @assert f([]) == 0    # f should be 0 at the empty set
  n = length(w)
  @assert length(S) == n
  i = sortperm(w, rev = true)
  V = sort(S)
  y = w[i]
  x = zeros(n)
  for ii = 1:n
    x[ii] = f(V[i[1: ii]]) - f(V[i[1: ii-1]])
  end
  return dot(x, y)
end
