#############################################################################
# frank_wolfe_away.jl
# Use the frank wolfe algorithm with away steps to solve problems of the DualModel:
#
# minimize g(x) st x ∈ P(F),
#
# where g(x) is a smooth convex function
# and P(F) is a polyhedron associated with a submodular function
#
# TODO:
# * improve data structure for ActiveVertices: would be faster to insert, delete, search with a heap
# * check maximization and minimization both work
#############################################################################

# import JuMP: DiffableFunction
import Base: deleteat!
export frank_wolfe_away

function frank_wolfe_away(p::SCOPEProblem{AssocPolyConstrained};
                          maxiters::Int64 = 100,
													epsilon::Float64 = 1e-3,
                          verbose::Bool = true)

  # initialize
  poly = p.model.poly
  n = length(poly.V)
  x = greedy(poly.F, collect(1:n))                  # the innitial solution
	g = grad(p.objective)
  V = ActiveVertices(length(x)) # maintain atomic representation in terms of vertices of polytope
	if verbose
		@printf("%10s%12s%10s\n", "iter", "obj", "\# active")
		@printf("%10d%12.4e%10d\n", 0, evaluate(p.objective, x), length(V.α))
	end

  for k = 1:maxiters

    gradient = evaluate(g, x)
    # compute fw step
    h = fenchel(poly, gradient)
		d_fw = h - x
		fw_dec = dot(-gradient, d_fw)     # the frank wolfe decrement
		# compute away step
		i, v, alpha = argmax_inner_product(gradient, V)
		d_away = x - v
		away_dec = dot(-gradient, d_away) # the away step decrement

		# check stopping condition
		if fw_dec < epsilon break end

    # choose away step or fw step
		if fw_dec > away_dec              # choose fw step
			d = d_fw
			gamma_max = 1
		else                              # choose away step
			d = d_away
			gamma_max = alpha / (1-alpha)
		end

    # choose stepsize gamma via backtracking linesearch with parameters (.1, .5)
		f(gamma) = evaluate(p.objective, x + gamma*d)
		# df(gamma) = dot(evaluate(g, x + gamma*d), d)
		df = dot(evaluate(g, x), d)

		gamma = gamma_max
		while f(gamma) > f(0) + .1 * df * gamma
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
				V = ActiveVertices([h], [1])
			else
				increase_weight!(V, h, gamma)
			end
		end
		if verbose && k%1==0
			@printf("%10d%12.4e%10d\n", 0, evaluate(p.objective, x), length(V.α))
		end
  end # maxiters

  return x
end # function

type ActiveVertices
	S::Array{Vector{Float64},1}
	α::Array{Float64,1}
end
ActiveVertices(n::Int) = ActiveVertices([zeros(n)], [1])

function argmax_inner_product(g::Vector, V::ActiveVertices)
	i = indmax([dot(g,s) for s in V.S])
	return i, V.S[i], V.α[i]
end

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
