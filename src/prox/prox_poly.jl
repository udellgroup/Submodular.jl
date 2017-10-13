#############################################################################
# prox_poly.jl
# Adapts the minimum-norm-point algorithm to compute the point in the submodular
# polyhedron p that is to an arbitrary given point w.
# For reference, look at "A Submodular Function Minimization Algorithm Based
# on the Minimum-Norm Base" by Fujishige and Isotani.
#############################################################################

export prox

TOL = 1e-3

function prox(p::AssocPoly, w::AbstractArray)
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
  projected = major_cycle(p, w, x, S)
  return projected[1]
end

# prox with warmstart
function prox(p::AssocPoly, w, x₀::AbstractArray)
  n = length(w)
  @assert length(p.V) == n
  S = zeros(n, 1)
  S[:, 1] = x₀
  projected = major_cycle(p, w, x₀, S)
  return projected[1]
end

function major_cycle(p::AssocPoly, w, x::AbstractArray, S::AbstractMatrix)
  n = length(w)
  S1 = minimum_check(p, w, x, S)
	if S1 == S
    return x
  else
    S₀ = [S; ones(1, size(S, 2))]
    x₀ = copy(x)
    push!(x₀, 1)
    λ = S₀\x₀                           # the coefficients of x in terms of points in S
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
        y = affproj(w, S)
        major_cycle(p, w, y, S)
      end
    else
      y = affproj(w, S)   # the point found
      y0 = copy(y)
      push!(y0, 1)
      S1 = [S; ones(1, size(S, 2))]
      μ = S1\y0                         # the coefficients of y in terms of points in S
      if all(μ .>= -TOL)
        major_cycle(p, w, y, S)
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
        major_cycle(p, w, x, S)
      end
    end
  end
end



function major_cycle(p::SubmodPoly, w, x::AbstractArray, S::AbstractMatrix)
  n = length(w)
  a = x - w
  if any(a .> 0)
    posind = find(a .> 0)
    # choose a new point to start iterations
    x₁ = x
    x₁[posind] = w[posind]
    z = fenchel(p, x₁ - w)
    b = z - w
    posind1 = find(b .> 0)
    z[posind1] = w[posind1]
    S1 = zeros(n, 1)
    S1[:, 1]  = z
    major_cycle(p, w, z, S1)
  else
    S1 = minimum_check(p, w, x, S)
		if S1 == S
	    return x
    else
			S₀ = [S; ones(1, size(S, 2))]
			x₀ = copy(x)
			push!(x₀, 1)
			λ = S₀\x₀                         # the coefficients of x in terms of points in S
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
					y = affproj(w, S)
					major_cycle(p, w, y, S)
				end
			else
				y = affproj(w, S)   # the point found
				y0 = copy(y)
				push!(y0, 1)
				S1 = [S; ones(1, size(S, 2))]
				μ = S1\y0                         # the coefficients of y in terms of points in S
				if all(μ .>= -TOL)
					major_cycle(p, w, y, S)
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
					major_cycle(p, w, x, S)
				end
			end
		end
	end
end

function minimum_check(p::AssocPoly, w, x::AbstractArray, S::AbstractMatrix)
  a = x - w
  z = fenchel(p, a)
  minprod = dot(z - w, a)
  curprod = dot(a, a)
  if curprod - minprod <= TOL
    return S
  else
    return [S z]
  end
end
