#############################################################################
# frank_wolfe_away.jl
# Use the frank wolfe algorithm with away steps to solve convex optimization of the form:
#
# minimize g(x) st x ∈ P(F),
#
# where g(x) is a smooth convex function
# and P(F) is a polyhedron associated with a submodular function
#
# todo:
# * improve data structure for ActiveVertices: would be faster to insert, delete, search with a heap
# * check maximization and minimization both work
#############################################################################

import JuMP: DiffableFunction
export frank_wolfe_away

function frank_wolfe_away(p::CombiProblem{DiffableFunction},
  maxiters=100,
  verbose=true)

  # check only constraint is one combinatorial constraint
  @assert length(p.convex_constraints==0) && length(p.combi_constraints==1)

  # initialize
  combiset = p.combi_constraints[1]
  n = combiset.dim
  x = zeros(n)
  V = ActiveVertices(length(x)) # maintain atomic representation in terms of vertices of polytope
	if verbose
		@printf("%10s%12s%10s\n", "iter", "obj", "\# active")
		@printf("%10d%12.4e%10d\n", 0, evaluate(p.objective, x), length(V.α))
	end

  for k=1:maxiters

    g = grad(p.objective, x)
    # compute fw step
    h = fenchel(combiset, g)
		d_fw = h - x
		fw_dec = dot(-g, d_fw) # the frank wolfe decrement
		# compute away step
		i, v, alpha = argmax_inner_product(-g, V)
		d_away = v - x
		away_dec = dot(-g, d_away) # the away step decrement

		# check stopping condition
		if fw_dec < epsilon break end

    # choose away step or fw step
		if fw_dec > away_dec
			d = d_fw
			gamma_max = 1
		else
			d = d_away
			gamma_max = alpha / (1-alpha)
		end

    # choose stepsize gamma via backtracking linesearch with parameters (.1, .5)
		f(gamma) = evaluate(p.objective, x + gamma*d)
		df(gamma) = dot(grad(p.objective, x + gamma*d), d)
		gamma = gamma_max
		while f(gamma) > f(0) + .1*df(0)*gamma
			gamma *= .5
		end

    # take the step and update atomic representation
		x += gamma*d
		if fw_dec <= away_dec # away step
			if verbose println("away! $(gamma/gamma_max)") end
			if gamma == gamma_max # drop step
				deleteat!(V, i)
			else # just decrease weight
				decrease_weight!(V, i, gamma)
			end
		else # FW step
			if verbose println("forward! $((gamma_max - gamma)/gamma_max)") end
			if gamma == 1 # we're at a vertex!
				V = ActiveVertices([s], [1])
			else
				increase_weight!(V, s, gamma)
			end
		end
		if verbose && k%1==0
			@printf("%10d%12.4e%10d\n", 0, evaluate(p.objective, x), length(V.α))
		end
  end # maxiters

  return x
end # function

type ActiveVertices
	S::Array{Vector{Int},1}
	α::Array{Float64,1}
end
ActiveVertices(n::Int) = ActiveVertices([zeros(Int, n)], [1])

function add!(V::ActiveVertices, α, v)
		push!(V.S, v)
		push!(V.α, α)
		z = sum(V.α) # normalization constant
		scale!(V.α, 1/z)
		return V
end
function deleteat!(V::ActiveVertices, i)
		deleteat!(V.S, i)
		deleteat!(V.α, i)
		z = sum(V.α) # normalization constant
		scale!(V.α, 1/z)
		return V
end
function decrease_weight!(V::ActiveVertices, i, gamma)
		scale!(V.α, 1+gamma)
		V.α[i] -= gamma
		return V
end
function increase_weight!(V::ActiveVertices, v, gamma)
		scale!(V.α, 1-gamma)
		idxs = find(a-> a==v, V.S)
		@assert length(idxs) <= 1 # there should never be replicates of the same vertex
		if length(idxs) > 0 # v is already in the set of active vertices
			V.α[idxs[1]] += gamma
		else # v is not yet in the set of active vertices
			push!(V.S, v)
			push!(V.α, gamma)
		end
		return V
end
function argmax_inner_product(g::Vector, V::ActiveVertices)
	i = indmax([dot(g,s) for s in V.S])
	return i, V.S[i], V.α[i]
end
