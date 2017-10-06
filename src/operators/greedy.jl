#############################################################################
# Computes the bases that maxize a linear function over the base polyhedron
# using the greedy algorithm.
#############################################################################

export greedy, greedy_rand

TOL = 1e-2

# Computes one solution
function greedy(f::Function, w::AbstractArray,
                S = collect(1:length(w))::AbstractArray)
  @assert f([])[1] == 0 # f should be 0 at the empty set
  n = length(w)
  @assert length(S) == n
  i = sortperm(w, rev = true)
  V = sort(S)
  x = zeros(n)
  for ii = 1:n
    x[i[ii]] = f(V[i[1: ii]])[1] - f(V[i[1: ii-1]])[1]
  end
  return x
end

greedy(f::CombiFunc, w::AbstractArray) = greedy(x -> evaluate(f, x), w)

# Compute a solution with indexed on the same level set permutated
function greedy_rand(f::Function, w::AbstractArray, tol = 1e-3 ::Number,
                    S = collect(1:length(w))::AbstractArray)
  @assert f([])[1] == 0       # f should be 0 at the empty set
  n = length(w)
  @assert length(S) == n
  ordering = sortperm(w, rev = true)
  pernum = 1
  V = sort(S)
  ordering = V[ordering]      # the indexed reordered in a descending order
  # Mark level sets
  index = 1
  level_sets = zeros(2, n)    # the beginning and ending of each level set
  level_sets[1, 1] = 1
  for i = 2:n
    if w[ordering[i - 1]] - w[ordering[i]] >= tol
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
  x = zeros(n)
  for ii = 1:n
    x[new_ordering[ii]] = f(V[new_ordering[1: ii]])[1] - f(V[new_ordering[1: ii-1]])[1]
  end
  return x, pernum, new_ordering
end

greedy_rand(f::CombiFunc, w::AbstractArray, tol::Number) = greedy_rand(x -> evaluate(f, x), w, tol)
