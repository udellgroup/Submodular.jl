#############################################################################
# Computes the bases that maxize a linear function over the base polyhedron
# using the greedy algorithm.
#############################################################################

export greedy, greedy_rand

TOL = 1e-2

# Computes one solution
function greedy(F::Function, x::AbstractArray,
                S = collect(1:length(x))::AbstractArray)
  @assert F([])[1] == 0 # f should be 0 at the empty set
  n = length(x)
  @assert length(S) == n
  i = sortperm(x, rev = true)
  V = sort(S)
  w = zeros(n)
  for ii = 1:n
    w[i[ii]] = F(V[i[1: ii]])[1] - F(V[i[1: ii-1]])[1]
  end
  return w
end

greedy(F::SubmodFunc, x::AbstractArray) = greedy(x -> evaluate(F, x), x)

greedy(f::LovaszExtAtom, x::AbstractArray) = greedy(f.func, x)

function greedy(f::LovaszExtAbsAtom, x::AbstractArray)
  w = greedy(f.func, abs.(x))
  for i = 1:length(w)
    if x[i] < 0
      w[i] = -w[i]
    end
  end
  return w
end

# Compute a solution with indexed on the same level set permutated
function greedy_rand(F::Function, x::AbstractArray, tol = 1e-3 ::Number,
                    S = collect(1:length(x))::AbstractArray)
  @assert F([])[1] == 0       # f should be 0 at the empty set
  n = length(x)
  @assert length(S) == n
  ordering = sortperm(x, rev = true)
  pernum = 1
  V = sort(S)
  ordering = V[ordering]      # the indexed reordered in a descending order
  # Mark level sets
  index = 1
  level_sets = zeros(2, n)    # the beginning and ending of each level set
  level_sets[1, 1] = 1
  for i = 2:n
    if x[ordering[i - 1]] - x[ordering[i]] >= tol
    level_sets[2, index] = i-1
    index += 1
    level_sets[1, index] = i
    end
  end
  level_sets[2, index] = n
  new_ordering = zeros(n)
  for i = 1:index
    new_ordering[Int(level_sets[1, i]) : Int(level_sets[2, i])] = shuffle(ordering[Int(level_sets[1, i]) : Int(level_sets[2, i])])
    pernum *= factorial(level_sets[2, i] - level_sets[1, i] + 1)
  end
  new_ordering = Int.(new_ordering)
  w = zeros(n)
  for ii = 1:n
    w[new_ordering[ii]] = F(V[new_ordering[1: ii]])[1] - F(V[new_ordering[1: ii-1]])[1]
  end
  return w, pernum, new_ordering
end

greedy_rand(F::SubmodFunc, x::AbstractArray, tol::Number) = greedy_rand(x -> evaluate(F, x), x, tol)
