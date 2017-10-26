#############################################################################
# minimum_norm_point.jl
# Implements the the minimum_norm_point algorithm, to find the point with
# minimum norm in the base polyhedron.
#############################################################################

export minimum_norm_point, prox

Tol = 1e-3

function minimum_norm_point(p::AssocPoly, w::AbstractArray)
  # step 1 initialization
  n = length(w)
  @assert length(p.V) == n
  V = sort(p.V)         # the sorted base set
  x = zeros(n)          # the candidate
  y = zeros(n + 1)
  for i = 1:n
    y[i + 1] = evaluate(p.f, V[1:i])[1]
    x[i] = y[i + 1] - y[i]
  end
  S = zeros(n, 1)       # the set of selected extreme points that constitute the affine set
  S[:, 1] = x
  # major cycle
  major_cycle(p, w, x, S)
end

function minimum_norm_point(p::PosPoly, w::AbstractArray)
  w₁ = copy(w)
  if any(w₁ .< 0)
    negind = find(w₁ .< 0)
    w₁[negind] = zeros(length(negind))
  end
  p₁ = SubmodPoly(p.f)
  minimum_norm_point(p₁, w₁)
end

function minimum_norm_point(p::SymPoly, w::AbstractArray)
  w₁ = sign.(w) .* w
  p₁ = SubmodPoly(p.f)
  x = minimum_norm_point(p₁, w₁)
  return sign.(w) .* x
end

# minimum-norm-point with warmstart
function minimum_norm_point(p::AssocPoly, w::AbstractArray, x₀::AbstractArray)
  n = length(w)
  @assert length(p.V) == n
  S = zeros(n, 1)
  S[:, 1] = x₀
  major_cycle(p, w, x₀, S)
end

function minimum_norm_point(p::PosPoly, w::AbstractArray, x₀::AbstractArray)
  w₁ = copy(w)
  if any(w₁ .< 0)
    negind = find(w₁ .< 0)
    w₁[negind] = zeros(length(negind))
  end
  p₁ = SubmodPoly(p.f)
  minimum_norm_point(p₁, w₁, x₀)
end

function minimum_norm_point(p::SymPoly, w::AbstractArray, x₀::AbstractArray)
  w₁ = sign.(w) .* w
  p₁ = SubmodPoly(p.f)
  x = minimum_norm_point(p₁, w₁, x₀)
  return sign.(w) .* x
end

function major_cycle(p::AssocPoly, w::AbstractArray, x::AbstractArray, S::AbstractMatrix)
  n = length(w)
  S1 = minimum_check(p, w, x, S)
	if S1 == S
    return x
  else
    S = S1
    colnum = size(S, 2)
    # step 3 find the point closes to w in the affine hull of points in S
    if colnum > size(S, 1)              # the affine hull is in face the whole space
      y0 = copy(w)
      push!(y0, 1)
      S1 = [S; ones(1, size(S, 2))]
      μ = S1\y0
      if all(μ .>= -Tol)                # μ is in the affine hull of points in set
        return w
      else                              # delete redundant points and compute the new x
        redind0 = find(μ .< -Tol)
        remind0 = deleteat!(collect(1:colnum), redind0)
        S = S[:, remind0]
        colnum = length(remind0)
        y = affproj(w, S)
        major_cycle(p, w, y, S)
      end
    else
      y = affproj(w, S)                 # the point found)
      y0 = copy(y)
      push!(y0, 1)
      S₁ = [S; ones(1, size(S, 2))]
      μ = S₁\y0                         # the coefficients of y in terms of points in S
      # if all(μ .>= -Tol)
      #   major_cycle(p, w, y, S)
      while any(μ .< -Tol)
        # step 4 find the closest point to y on the segment between x and y while staying in the convex hull of points in S, and set it to be x
        x₀ = copy(x)
        push!(x₀, 1)
        λ = S₁\x₀                           # the coefficients of x in terms of points in S
        ratios = λ ./ μ
        negratios = ratios[find(μ .< -Tol)]
        α = maximum(negratios)
        redind = find((ratios .<= α + Tol) .& (ratios .<= α + Tol))    # the index of the redundant points
        remind = deleteat!(collect(1:colnum), redind)
        S = S[:, remind]                           # the new set of points
        x = y                                      # the new candidate
        y = affproj(w, S)                          # the new optimal solution given S
        y0 = copy(y)
        push!(y0, 1)
        S₁ = [S; ones(1, size(S, 2))]
        μ = S₁\y0
      end
      major_cycle(p, w, y, S)
    end
  end
end

function minimum_check(p::AssocPoly, w::AbstractArray, x::AbstractArray, S::AbstractMatrix)
  a = x - w
  z = fenchel(p, a)
  minprod = dot(z - w, a)
  curprod = dot(a, a)
  if curprod - minprod <= Tol
    return S
  else
    return [S z]
  end
end

function minimum_check(p::SubmodPoly, w::AbstractArray, x::AbstractArray, S::AbstractMatrix)
  a = x - w
  pos = false
  if any(a .> 0)
    pos = true
    posind = find(a .> 0)
  end
  z = greedy(p.f, -a)
  if pos                                                 # x - w has positive elements
    s = z[posind] - w[posind]
    if any(s .> 0)                                       # some elements of z should be adjusted
      posind1 = find(s .> 0)
      z[posind[posind1]] = w[posind[posind1]]
    end
  end
  minprod = dot(z - w, a)
  curprod = dot(a, a)
  if curprod - minprod <= Tol
    return S
  else
    S1 = [S z]
    return [S z]
  end
end

function minimum_norm_point(f::SubmodFunc, w::AbstractArray)
  p = BasePoly(f)
  return minimum_norm_point(p, w)
end

prox(p::AssocPoly, w::AbstractArray) = minimum_norm_point(p, w)
