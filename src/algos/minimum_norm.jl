#############################################################################
# minimum_norm.jl
# Adapts the minimum-norm-point algorithm to compute the point in the base
# polyhedron p that is to an arbitrary given point w.
# For reference, look at "A Submodular Function Minimization Algorithm Based
# on the Minimum-Norm Base" by Fujishige and Isotani.
#############################################################################

export minimum_norm

TOL = 1e-6

function minimum_norm(f::Function, T::AbstractArray, w::AbstractVector)
  # step 1 initialization
  n = length(w)
  @assert length(T) == n
  V = sort(T)       # the sorted base set
  x = zeros(n)      # the candidate
  for i = 1:n
    x[i] = f(V[1:i]) - f(V[1:i-1])
  end
  S = zeros(n, 1)   # the set of selected extreme points that constitute the affine set
  S[:, 1] = x
  # major cycle
  major_cycle(f, T, w, x, S)
end

# step 2 find the point z that minimize (x-w)^(z-w) st z in the base polyhedron
function minimum_check(f::Function, V::AbstractArray, w, x::AbstractVector, S::AbstractMatrix)
  a = x - w
  z = greedy(f, -a, V)
  minprod = dot(z - w, a)
  curprod = dot(a, a)
  if curprod - minprod < TOL
    return S
  else
    return [S z]
  end
end

function major_cycle(f::Function, V::AbstractArray, w, x::AbstractVector, S::AbstractMatrix)
  n = length(w)
  S1 = minimum_check(f, V, w, x, S)
  if S1 == S
    return x
  else
    S0 = [S; ones(1, size(S, 2))]
    x0 = copy(x)
    push!(x0, 1)
    λ = S0\x0                         # the coefficients of x in terms of points in S
    push!(λ, 0)
    S = S1
    colnum = size(S, 2)
    # step 3 find the point closes to w in the affine hull of points in S
    if colnum > size(S, 1)              # the affine hull is in face the whole space
      y0 = copy(w)
      push!(y0, 1)
      S1 = [S; ones(1, size(S, 2))]
      μ = S1\y0
      if all(μ .>= -TOL)                # μ is in the affine hull of points in set
        return w
      else                              # delete redundant points and compute the new x
        redind0 = find(μ .< -TOL)
        remind0 = deleteat!(collect(1:colnum), redind0)
        S = S[:, remind0]
        colnum = length(remind0)
        y = affpro(w, S)
        major_cycle(f, V, w, y, S)
      end
    else
      y = affpro(w, S)   # the point found
      y0 = copy(y)
      push!(y0, 1)
      S1 = [S; ones(1, size(S, 2))]
      μ = S1\y0                         # the coefficients of y in terms of points in S
      if all(μ .>= -TOL)
        major_cycle(f, V, w, y, S)
      else
        # step 4 find the closest point to y on the segment between x and y while
        # staying in the convex hull of points in S, and set it to be x
        ratios = λ ./ μ
        negratios = ratios[find(μ .< -TOL)]
        α = maximum(negratios)
        redind = find((ratios .<= α + TOL) .& (ratios .<= α + TOL))    # the index of the redundant points
        remind = deleteat!(collect(1:colnum), redind)
        S = S[:, remind]                           # the new set of points
        β = α/(α - 1)
        x = β * y + (1 - β) * x                   # the new candidate
        major_cycle(f, V, w, x, S)
      end
    end
  end
end

function minimum_norm(f::Function, n::Int, w::AbstractVector)
  return minimum_norm(f::Function, collect(1:n), w::AbstractVector)
end
