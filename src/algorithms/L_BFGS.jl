#############################################################################
# L_BFGS.jl
# Implement the L_BFGS method to solve convex optimization of the form:
# minimize f(x) + g(x), where f(x) is the Lovasz extension of a submodular funciton.
#############################################################################

export L_BFGS

function L_BFGS(p::Problem,                    # the convex problem
                f::LovaszExtAtom;              # the lovasz extension
                x₀::AbstractArray = collect(1:length(get_v(f)[1])),  # the starting point
                m::Int = length(get_v(f)[1]),  # the number of kept iterations
                c₁::Float64 = 0.7,             # linesearch parameter
                c₂::Float64 = 0.9,             # linesearch parameter
                μ::Float64 = 0.9,              # linesearch shrink rate
                epsilon::Float64 = 1e-3,       # the stopping criteria of the algorithm
                max_iters::Int = 100)          # the maximal iteration time

  if length(get_v(p.objective + f)) != 1
    error("Cannot optimize functions with multiple variables.")
  end

  x = x₀                                              # the candidate
  objective = p.objective + f
  var = get_v(f)[1]
  n = length(var)
  s = Array{AbstractArray}(m)                         # the differences between candidates of neighbouring iterations
  y = Array{AbstractArray}(m)                         # the differences between gradients candidates of neighbouring iterations
  α = Array{Number}(m)
  ρ = Array{Number}(m)
  r = AbstractArray                                   # quasi_gradient
  z = AbstractArray                                   # descent direction
  λ = 1.0                                             # step size
  λ̃ = 1.0                                             # the λ producing the largest gap satisfying the wolfe condition
  β = 0
  decr = 0.0                                          # the decrement
  gap = 0.0
  gap_upper = 0.0
  check = false
  grad_conv = grad(p.objective)
  obj(λ) = evaluate(objective, x + λ*z)[1]
  gra(λ) = evaluate(grad_conv, x + λ*z) + greedy(f.func, x + λ*z)
  q = evaluate(grad_conv, x) + greedy(f.func, x)

  for k = 0:max_iters
    # println("k = $k")
    # println("norm(q) = $(norm(q))")
    if norm(q) < epsilon
      return x
      break
    else
      # lev_tol = norm(q)/(n^0.5)*0.01                  # the tolerance for level sets
      if k == 0
        r = q
      elseif k < m
        H = (dot(s[k], y[k]) / dot(y[k], y[k])) * eye(n)
        for i = k:-1:1
          α[i] = ρ[i] * dot(s[i], q)
          q -= α[i] * y[i]
        end
        r = H * q
        for i = 1:k
          β = ρ[i] * dot(y[i], r)
          r += s[i] * (α[i] - β)
        end
      else
        k̃ = mod(k-1, m) + 1
        # println("k̃ = $k̃")
        H = (dot(s[k̃], y[k̃]) / dot(y[k̃], y[k̃])) * eye(n)
        for i = k:-1:k-m+1
          ĩ = mod(i-1, m) + 1
          α[ĩ] = ρ[ĩ] * dot(s[ĩ], q)
          q -= α[ĩ] * y[ĩ]
        end
        r = H * q
        for i = (k-m+1):k
          ĩ = mod(i-1, m) + 1
          β = ρ[ĩ] * dot(y[ĩ], r)
          r += s[ĩ] * (α[ĩ] - β)
        end
      end
      z = -r
      λ = 1
      q = evaluate(grad_conv, x) + greedy(f.func, x)
      obj0 = obj(0)
      decr = dot(q, z)
      # println("x = $x")
      # println("q = $q")
      # println("z = $z")
      # println("decr = $decr")
      if abs(decr) < epsilon
        return x
        break
      end
      while obj(λ) > (obj0 + c₁ * λ * decr)
        if λ < 1e-2
          error("Line search fail to terminate")
          return x
        end
        # println("λ = $λ")
        # println("cond1 = $(obj(λ) - (obj0 + c₁ * λ * decr))")
  		  λ *= μ
  	  end
      while dot(gra(λ), z) < c₂ * decr
        if λ < 1e-2
          error("Line search fail to terminate")
          return x
        end
        # println("λ = $λ")
        # println("cond2 = $(dot(gra(λ), z) - c₂ * decr)")
  		  λ *= μ
      end
      # λ̃ = copy(λ)
      # check, gap = level_detector(x + λ*z, lev_tol)
      # gap_upper = gap
      # while check == true && dot(gra(λ), z) >= c₂ * decr
      #   if gap_upper < gap
      #     λ̃ = λ
      #     gap_upper = gap
      #   end
      #   λ *= μ
      #   # println(λ)
      #   check, gap = level_detector(x + λ*z, lev_tol)
      # end
      # if dot(gra(λ), z) < c₂ * decr                   # the wolfe condition is violated, pick the λ that produces the largest gap while satisfying the wolfe condition
      #   λ = λ̃
      # end
      # println()
      # println("λ = $λ")
      # println()
      x += λ * z
      k̃ = mod(k, m) + 1                              # update new entries
      s[k̃] = λ * z
      q_new = evaluate(grad_conv, x) + greedy(f.func, x)
      y[k̃] = q_new - q
      ρ[k̃] = 1 / dot(y[k̃], s[k̃])
      q = q_new
    end
  end
end

function level_detector(x::AbstractArray, lev_tol::Float64 = 1e-2)
  w = sort(x)
  gap = lev_tol                                       # the tightest gap
  check = false
  for i = 1:(length(x)-1)
    if w[i+1] - w[i] < gap
      gap = w[i+1] - w[i]
      check = true
    end
  end
  return check, gap
end
