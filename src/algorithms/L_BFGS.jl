#############################################################################
# L_BFGS.jl
# Implement the L_BFGS method to solve convex optimization of the form:
# minimize f(x) + g(x), where f(x) is the Lovasz extension of a submodular funciton.
#############################################################################

export L_BFGS

function L_BFGS(p::Problem,                    # the convex problem
                f::LovaszExtAtom;              # the lovasz extension
                x₀::AbstractArray = zeros(length(get_v(f)[1])),  # the starting point
                m::Int = length(get_v(f)[1]),  # the number of kept iterations
                c₁::Float64 = 0.5,             # linesearch parameter
                c₂::Float64 = 0.1,             # linesearch parameter
                μ::Float64 = 0.7,              # linesearch shrink rate
                epsilon::Float64 = 1e-3,       # the stopping criteria of the algorithm
                max_iters::Int = 100)          # the maximal iteration time

  if length(get_v(p.objective + f)) != 1
    error("Cannot optimize functions with multiple variables.")
  end

  x = x₀                                              # the candidate
  var = get_v(f)[1]
  n = length(var)
  s = Array{AbstractArray}(m)                         # the differences between candidates of neighbouring iterations
  y = Array{AbstractArray}(m)                         # the differences between gradients candidates of neighbouring iterations
  α = Array{Number}(m)
  ρ = Array{Number}(m)
  r = AbstractArray                                   # descent direction
  α = 1                                               # step size
  β = 0
  g_convex = grad(p.objective)
  q = evaluate(g_convex, x) + greedy(f.func, x)

  for k = 0:max_iters
    println(k)
    if norm(q) < epsilon
      return x
      break
    else
      if k == 0
        r = q
      else
        H = (dot(s[k], y[k]) / dot(y[k], y[k])) * eye(n)
        if k < m
          for i = k:-1:1
            println("α = $(α[i])")
            println("ρ = $(ρ[i])")
            println("s = $(s[i])")
            println("q = $q")
            α[i] = ρ[i] * dot(s[i], q)
            println("53")
            q -= α[i] * y[i]
          end
          r = H * q
          for i = 1:k
            β = ρ * y[i] * r
            r += s[i] * (α[i] - β)
          end
        else
          for i = k:-1:k-m+1
            ĩ = mod(i-1, m) + 1
            α[ĩ] = ρ[ĩ] * dot(s[ĩ], q)
            q -= α[ĩ] * y[ĩ]
          end
          r = H * q
          for i = (k-m+1):k
            ĩ = mod(i-1, m) + 1
            β = ρ * y[ĩ] * r
            r += s[ĩ] * (α[ĩ] - β)
          end
        end
      end
      z = -r
      α = 1
      q = evaluate(g_convex, x) + greedy(f.func, x)
      obj(α) = evaluate(p.objective + f, x + α*z)[1]
      gra(α) = evaluate(g_convex, x + α*z)[1] + greedy(f.func, x + α*z)
      while obj(α) > obj(0) + c₁ * α * dot(q, z) || dot(gra(α), z) < c₂ * dot(q, z)
  		  α *= μ
  	  end
      x += α * z
      k̃ = mod(k, m) + 1                              # update new entries
      s[k̃] = α * z
      q_new = evaluate(g_convex, x) + greedy(f.func, x)
      y[k̃] = q_new - q
      ρ[k̃] = 1 / dot(y[k̃], s[k̃])
      q = q_new
    end
  end
end
