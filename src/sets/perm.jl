import Base.in
export in, Perm, fenchel


# the set of permutations of 1,...,n
type Perm #<:CombiSet
  dim::Int
end

# check if x is in the Permutation set
function in(c::AbstractVector, p::Perm)
  n = p.dim
  @assert length(c) == n
  return Set(c) == Set(collect(1:n))
end

# compute min c^x st x in p
function fenchel(p::Perm, c::AbstractVector)
  n = p.dim
  @assert length(c) == n
  i = sortperm(c, rev = true) # biggest index of c is first element of i, etc.
  h = zeros(n)
  for ii=1:n
    h[i[ii]] = ii
  end
  return h
end
