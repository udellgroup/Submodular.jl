#############################################################################
# proximal_level_bundle.jl
# Implements the proximal level bundle method.
# For reference look at:
#############################################################################

export proximal_level_bundle

TOL = 1e-3

# maxiters = 20

# λ = 0.5 * 2^0.5
λ = 1

function proximal_level_bundle(p::Problem, f::LovaszExtAtom)
  if length(get_v(p.objective)) > 1
    error("Cannot optimize functions with multiple variables.")
  end

  # the approximated problem
  q = Problem(p.head, p.objective, p.constraints)
  obj = p.objective + f                          # the objective function

  # warm start
  solve!(q, SCSSolver(verbose=false))
  t = Variable(1)
  q.objective += t
  push!(q.solution.primal, q.solution.primal[end])
  q.solution.primal[end - 1] = 0
  q.solution.status = "not yet solved"
  var = get_v(q.objective.children[1])[1]
  a = greedy_rand(f.func, var.value[:])
  push!(q.constraints, dot(a, var) <= t)
  approx = []                                    # the affine functions
  push!(approx, q.constraints[end].lhs)

  # Step 0: Initialization
  solve!(q, SCSSolver(verbose=false), warmstart = true)
  t.value = evaluate(approx[1])[1]
  upper = evaluate(obj)[1]
  lower = evaluate(q.objective)[1]
  gap = upper - lower
  a = greedy_rand(f.func, var.value[:])
  push!(q.constraints, dot(a, var) <= t)
  push!(approx, q.constraints[end].lhs)

  # level, set as a variable so its value updates automatically
  level = Variable(1)
  fix!(level, λ * lower + (1 - λ) * upper)       # the level set

  # projection center
  size = var.size
  x₀ = Variable(size)                            # the prox center
  fix!(x₀, var.value[:])

  # projection
  proj = Problem(:minimize, norm(var - x₀))
  dotprod = dot(a, var)
  diff = evaluate(obj)[1] - evaluate(dotprod)[1]
  push!(proj.constraints, dotprod + diff - level <= 0)
  pos_test = []
  push!(pos_test, proj.constraints[end].lhs)
  println("diff = $(evaluate(pos_test[1])[1])")

  iters = 0                                      # TODO: update TODO: k(l)
  approx_iters = 0
  optimal = var.value[:]                         # the current optimal solution
  cur_val = 0
  pos(x) = evaluate(x)[1] > -1e-1                # find positive elements

  for i = 1:maxiters
    println("gap = $gap")
    println()
    # Step 1: optimality check
    if gap < TOL
      println(iters)
      println(approx_iters)
      return optimal
      break                                      # optimality reached
    else
      # Step 2: level feasibility check
      solve!(proj, SCSSolver(verbose=false, suppress_warnings=true), warmstart = iters > 0)
      iters += 1
      println("proj.optval = $(proj.optval)")
      if proj.status == :Infeasible
        # Step 3: update lower bound and level
        println("lower bound update")
        solve!(q, SCSSolver(verbose=false), warmstart = true)
        approx_iters += 1
        t.value = minimum(evaluate.(approx))[1]
        lower = evaluate(q.objective)[1]
        level.value = λ * lower + (1 - λ) * upper
      elseif proj.status == :Unbounded
        break
        error("The convex function is not strongly convex")
      else
        # Step 4: update linear constraints
        # start_index = length(proj.solution.dual) - length(approx)
        # pos_index = find(pos_test, proj.solution.dual[(start_index + 1): end]) # update the trial set
        maxx = -10000
        for i = 1:length(pos_test)
          maxx = max(maxx, evaluate(pos_test[i])[1])
        end
        println("maxx = $maxx")
        pos_index = find(pos, pos_test)        # update the trial set
        println("pox_index = $pos_index")
        q.constraints = q.constraints[pos_index]
        approx = approx[pos_index]
        proj.constraints = proj.constraints[pos_index]
        pos_test = pos_test[pos_index]
        a = greedy_rand(f.func, var.value[:])
        push!(q.constraints, dot(a, var) <= t)
        push!(approx, q.constraints[end].lhs)
        dotprod = dot(a, var)
        diff = evaluate(obj)[1] - evaluate(dotprod)[1]
        push!(proj.constraints, dotprod + diff - level <= 0)
        push!(pos_test, proj.constraints[end].lhs)
        maxx1 = -10000
        for i = 1:length(pos_test)
          maxx1 = max(maxx1, evaluate(pos_test[i])[1])
        end
        println("maxx1 = $maxx1")
        # Step 5: add new candidate and update upper bound and level
        x₀.value = var.value[:]                # update the prox centers
        cur_val = evaluate(obj)[1]
        if upper > cur_val
          optimal = var.value[:]
          upper = cur_val
          gap = upper - lower
          level.value = λ * lower + (1 - λ) * upper
        # else
        #   var.value[:] = optimal
        end
      end
    end
  end
end
