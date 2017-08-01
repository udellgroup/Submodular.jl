# compute the base that maxizes a linear function over the base polyhedron using the greedy algorithm

export greedy

function greedy(f::Function, w::AbstractVector,
                S = collect(1:length(w))::AbstractVector)
  @assert f([]) == 0 # f should be 0 at the empty set
  n = length(w)
  @assert length(S) == n
  i = sortperm(w, rev = true)
  V = sort(S)
  x = zeros(n)
  for ii = 1:n
    x[ii] = f(V[i[1: ii]]) - f(V[i[1: ii-1]])
  end
  return x
end
