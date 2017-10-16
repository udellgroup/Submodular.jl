#############################################################################
# cutting_plane.jl
# Use the cutting plane method to solve convex optimization of the form:
# minimize f(x) + g(x), where f(x) is the Lovasz extension of a submodular funciton.
#############################################################################

# using SCS
# using ECOS
#using Gurobi
using Mosek

export cutting_plane

function cutting_plane(p::Problem, f::LovaszExtAtom, opt_tol, max_iters, max_rep, lev_tol_def, lev_tol_up, lev_tol_low, shrink_point, λ, μ)

  lev_tol = lev_tol_def

  if length(get_v(p.objective + f)) != 1
    error("Cannot optimize functions with multiple variables.")
  end
  q = Problem(p.head, p.objective, p.constraints)
  conslength = length(p.constraints)

  # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = 1e-3))
  # solve!(q, SCSSolver(verbose=false))
  # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
  solve!(q, MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-3 * opt_tol))

  # Reset the primal variable for warmstart
  push!(q.solution.primal, q.solution.primal[end])
  q.solution.primal[end - 1] = 0

  var = get_v(q.objective)[1]
  t = Variable(1)
  q.objective += t
  n = size(var)[1]
  a = greedy(f.func, var.value[:])
  push!(q.constraints, dot(a, var) <= t)
  constraints = []
  push!(constraints, q.constraints[end].lhs)       # the array of the affine functions
  subgradients = [a]

  cur_upper = 0.0
  cur_low = 0.0
  upper = Inf
  lower = -Inf
  optsol = var.value[:]
  pernum = 1
  permset = [sortperm(optsol, rev = true)]         # orderings explored
  low_sol = zeros(length(var.value[:]))
  iters = 1
  rep = 0                                          # the repetition count of a trial point
  cos_num = 0.0
  gap = 0.0

  lastt = Problem(q.head, q.objective, q.constraints)

  qq = Problem(q.head, q.objective, q.constraints)

  for i = 1:max_iters

    # println(length(permset))

    # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = TOL^2), warmstart = true)
    # solve!(q, SCSSolver(verbose=false, suppress_warnings = true), warmstart = true)
    # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
    solve!(q, MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-3 * opt_tol))

    if isnan(q.optval)
      var.value[:] = optsol
      # println("Mosek returns NaN!")
      # println("gapopt = $gap")
      break
    end

    lastt = Problem(q.head, q.objective, q.constraints)

    inacct = t.value                             # the inaccurate value of t
    acct = maximum(evaluate.(constraints))       # the accurate value of t
    cur_sol = var.value[:]
    cur_upper = q.optval - inacct + evaluate(f)[1]
    if cur_upper < upper
      optsol = var.value[:]
      upper = copy(cur_upper)
    end

    if q.status == :Optimal
      cur_low = q.optval - inacct + acct
    else
      cur_low = solvedualMosek!(q)
    end

    if lower < cur_low                           # update the lower bound
      low_sol = var.value[:]
    # else
      # println("lowstab")
    end
    lower = max(cur_low, lower)
    gap = upper - lower
    if gap < opt_tol
      var.value[:] = optsol
      # println("gapopt = $gap")
      break
    # else
      # println("gapsub = $gap")
    end
    if upper < cur_upper
      rep += 1
      if length(permset) >= max_rep
        var.value[:] = optsol
        println("reached the maximum iteration at a single point")
        # println("gapsub = $gap")
        break
      end
      if rep == shrink_point                         # check if the solver is stuck
        while pernum > max_rep && level_tol >= lev_tol_low  # println("reduce permnum!")
          lev_tol *= μ
          (c, pernum, perm) = greedy_rand(f.func, optsol, lev_tol)
        end
      end
      b = greedy(f.func, var.value[:])
      while in(b, subgradients)
        # println("b included!")
        (b, pernum1, perm) = greedy_rand(f.func, var.value[:], lev_tol)
      end
      var.value[:] = optsol
      while length(permset) >= pernum && level_tol <= lev_tol_up  # println("enlarge permnum!")
        lev_tol *= λ
        (a, pernum, perm) = greedy_rand(f.func, optsol, lev_tol)
      end
      (a, pernum, perm) = greedy_rand(f.func, optsol, lev_tol)
      while in(a, subgradients)
        # println("a included!")
        if !in(perm, permset)
          push!(permset, perm)
        end
        (a, pernum) = greedy_rand(f.func, optsol, lev_tol)
      end
      push!(permset, perm)
    else                                           # update the upper bound
      rep = 1                                      # reset repetition count
      lev_tol = lev_tol_def                        # reset the level set tolerance
      (a, pernum, perm) = greedy_rand(f.func, optsol, lev_tol)
      a = greedy(f.func, optsol)
      permset = [sortperm(var.value[:], rev = true)]
    end

    cos_num = length(constraints)
    if cos_num > n
      prod = zeros(cos_num)
      for i = 1:cos_num
        prod[i] = dot(optsol, subgradients[i])
      end
      if cos_num == (n + 2)
        minn = indmin(prod)
        minn = indmin(prod)
        q.constraints[conslength + minn] = q.constraints[end]
        constraints[minn] = constraints[end]
        subgradients[minn] = subgradients[end]
        prod = prod[1:end-1]
        q.constraints = q.constraints[1:end-1]
        constraints = constraints[1:end-1]
        subgradients = subgradients[1:end-1]
      end
      if rep == 1
        minn = indmin(prod)
        q.constraints[conslength + minn] = (dot(a, var) <= t)
        constraints[minn] = q.constraints[conslength + minn].lhs
        subgradients[minn] = a
      else                                       # b is added
        minn = indmin(prod)
        q.constraints[conslength + minn] = (dot(a, var) <= t)
        constraints[minn] = q.constraints[conslength + minn].lhs
        subgradients[minn] = a
        maxx = indmin(prod)
        prod[minn] = prod[maxx]
        minn1 = indmin(prod)
        q.constraints[conslength + minn1] = (dot(b, var) <= t)
        constraints[minn1] = q.constraints[conslength + minn1].lhs
        subgradients[minn1] = b
      end
    else
      if rep > 1
        push!(q.constraints, dot(b, var) <= t)
        push!(constraints, q.constraints[end].lhs)
        push!(subgradients, b)
      end
      push!(q.constraints, dot(a, var) <= t)
      push!(constraints, q.constraints[end].lhs)
      push!(subgradients, a)
    end
    iters += 1
  end
  return optsol
end

function cutting_plane(p::Problem, f::LovaszExtAbsAtom, opt_tol, max_iters, max_rep, lev_tol_def, lev_tol_up, lev_tol_low, shrink_point, λ, μ)

  lev_tol = lev_tol_def

  if length(get_v(p.objective + f)) != 1
    error("Cannot optimize functions with multiple variables.")
  end
  q = Problem(p.head, p.objective, p.constraints)
  conslength = length(p.constraints)

  # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = 1e-3))
  # solve!(q, SCSSolver(verbose=false))
  # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
  solve!(q, MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-3 * opt_tol))

  # Reset the primal variable for warmstart
  push!(q.solution.primal, q.solution.primal[end])
  q.solution.primal[end - 1] = 0

  var = get_v(q.objective)[1]
  t = Variable(1)
  q.objective += t
  n = size(var)[1]
  a = greedy(f.func, abs.(var.value[:]))
  push!(q.constraints, dot(a, abs(var)) <= t)
  constraints = []
  push!(constraints, q.constraints[end].lhs)       # the array of the affine functions
  subgradients = [a]

  cur_upper = 0.0
  cur_low = 0.0
  upper = Inf
  lower = -Inf
  optsol = var.value[:]
  pernum = 1
  permset = [sortperm(abs.(optsol), rev = true)]         # orderings explored
  low_sol = zeros(length(var.value[:]))
  iters = 1
  rep = 0                                          # the repetition count of a trial point
  cos_num = 0.0
  gap = 0.0

  lastt = Problem(q.head, q.objective, q.constraints)

  qq = Problem(q.head, q.objective, q.constraints)

  for i = 1:max_iters

    # println(length(permset))

    # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = TOL^2), warmstart = true)
    # solve!(q, SCSSolver(verbose=false, suppress_warnings = true), warmstart = true)
    # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
    solve!(q, MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-3 * opt_tol))

    if isnan(q.optval)
      var.value[:] = optsol
      # println("Mosek returns NaN!")
      # println("gapopt = $gap")
      break
    end

    lastt = Problem(q.head, q.objective, q.constraints)

    inacct = t.value                             # the inaccurate value of t
    acct = maximum(evaluate.(constraints))       # the accurate value of t
    cur_sol = var.value[:]
    cur_upper = q.optval - inacct + evaluate(f)[1]
    if cur_upper < upper
      optsol = var.value[:]
      upper = copy(cur_upper)
    end

    if q.status == :Optimal
      cur_low = q.optval - inacct + acct
    else
      cur_low = solvedualMosek!(q)
    end

    if lower < cur_low                           # update the lower bound
      low_sol = var.value[:]
    # else
      # println("lowstab")
    end
    lower = max(cur_low, lower)
    gap = upper - lower
    if gap < opt_tol
      var.value[:] = optsol
      # println("gapopt = $gap")
      break
    # else
      # println("gapsub = $gap")
    end
    if upper < cur_upper
      rep += 1
      if length(permset) >= max_rep
        var.value[:] = optsol
        println("reached the maximum iteration at a single point")  # println("gapsub = $gap")
        break
      end
      if rep == shrink_point                         # check if the solver is stuck
        while pernum > max_rep && level_tol >= lev_tol_low  # println("reduce permnum!")
          lev_tol *= μ
          (c, pernum, perm) = greedy_rand(f.func, abs.(optsol), lev_tol)
        end
      end
      b = greedy(f.func, abs.(var.value[:]))
      while in(b, subgradients)
        # println("b included!")
        (b, pernum1, perm) = greedy_rand(f.func, abs.(var.value[:]), lev_tol)
      end
      var.value[:] = optsol
      while length(permset) >= pernum && level_tol <= lev_tol_up  # println("enlarge permnum!")
        lev_tol *= λ
        (a, pernum, perm) = greedy_rand(f.func, abs.(optsol), lev_tol)
      end
      (a, pernum, perm) = greedy_rand(f.func, abs.(optsol), lev_tol)
      while in(a, subgradients)
        # println("a included!")
        if !in(perm, permset)
          push!(permset, perm)
        end
        (a, pernum) = greedy_rand(f.func, abs.(optsol), lev_tol)
      end
      push!(permset, perm)
    else                                           # update the upper bound
      rep = 1                                      # reset repetition count
      lev_tol = lev_tol_def                        # reset the level set tolerance
      (a, pernum, perm) = greedy_rand(f.func, abs.(optsol), lev_tol)
      a = greedy(f.func, abs.(optsol))
      permset = [sortperm(var.value[:], rev = true)]
    end

    cos_num = length(constraints)
    if cos_num > n
      prod = zeros(cos_num)
      for i = 1:cos_num
        prod[i] = dot(optsol, subgradients[i])
      end
      if cos_num == (n + 2)
        minn = indmin(prod)
        minn = indmin(prod)
        q.constraints[conslength + minn] = q.constraints[end]
        constraints[minn] = constraints[end]
        subgradients[minn] = subgradients[end]
        prod = prod[1:end-1]
        q.constraints = q.constraints[1:end-1]
        constraints = constraints[1:end-1]
        subgradients = subgradients[1:end-1]
      end
      if rep == 1
        minn = indmin(prod)
        q.constraints[conslength + minn] = (dot(a, abs(var)) <= t)
        constraints[minn] = q.constraints[conslength + minn].lhs
        subgradients[minn] = a
      else                                       # b is added
        minn = indmin(prod)
        q.constraints[conslength + minn] = (dot(a, abs(var)) <= t)
        constraints[minn] = q.constraints[conslength + minn].lhs
        subgradients[minn] = a
        maxx = indmin(prod)
        prod[minn] = prod[maxx]
        minn1 = indmin(prod)
        q.constraints[conslength + minn1] = (dot(b, abs(var)) <= t)
        constraints[minn1] = q.constraints[conslength + minn1].lhs
        subgradients[minn1] = b
      end
    else
      if rep > 1
        push!(q.constraints, dot(b, abs(var)) <= t)
        push!(constraints, q.constraints[end].lhs)
        push!(subgradients, b)
      end
      push!(q.constraints, dot(a, abs(var)) <= t)
      push!(constraints, q.constraints[end].lhs)
      push!(subgradients, a)
    end
    iters += 1
  end
  return optsol
end
