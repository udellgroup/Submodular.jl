# the set of permutations of 1,...,n
type Permutation #<:CombiSet
  dim::Int
end

# compute min c^x st x in p
function fenchel(p::Permutation, c::AbstractVector)
  n = p.dim
  @assert length(c) == n
  i = sortperm(c) # smallest index of c is first element of i, etc. this is the inverse of the permutation we want.
  h = zeros(n)
  for ii=1:n
    h[i[ii]] = n - ii + 1
  end
  return h
end
