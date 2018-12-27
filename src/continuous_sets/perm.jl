import Base.in
export in, Perm, fenchel


# the set of permutations of 1,...,n
mutable struct Perm <: ContiSet
  dim::Int
end

# check if w is in the Permutation set
function in(w::AbstractArray, p::Perm)
  n = p.dim
  @assert length(w) == n
  return Set(w) == Set(collect(1:n))
end

# compute min w^x st x in p
function fenchel(p::Perm, w::AbstractArray)
  n = p.dim
  @assert length(w) == n
  i = sortperm(w, rev = true) # biggest index of c is first element of i, etc.
  h = zeros(n)
  for ii=1:n
    h[i[ii]] = ii
  end
  return h
end
