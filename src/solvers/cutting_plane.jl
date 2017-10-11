#############################################################################
# cutting_plane.jl
# Use the cutting plane method to solve convex optimization of the form:
# minimize f(x) + g(x), where f(x) is the Lovasz extension of a submodular funciton.
#############################################################################

# using SCS
using ECOS
#using Gurobi
using Mosek

export cutting_plane, cutting_plane1

TOL1 = 1e-3

tight_tol = 1e-3

soft_tol = 1e-2

maxperm = 100

checkpoint = 20

maxiters = 200

function cutting_plane(p::Problem, f::LovaszExtAtom)

  lsTOL = 1e-3

  if length(get_v(p.objective + f)) != 1
    error("Cannot optimize functions with multiple variables.")
  end
  q = Problem(p.head, p.objective, p.constraints)
  # obj = p.objective + f
  conslength = length(p.constraints)

  # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = 1e-3))
  # solve!(q, SCSSolver(verbose=false))
  # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
  solve!(q, MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-6))

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
  haha = []

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
  max_rep = 0

  lastt = Problem(q.head, q.objective, q.constraints)

  qq = Problem(q.head, q.objective, q.constraints)

  for i = 1:maxiters

    println(length(permset))

    # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = TOL^2), warmstart = true)
    # solve!(q, SCSSolver(verbose=false, suppress_warnings = true), warmstart = true)
    # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
    solve!(q, MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-6))

    if isnan(q.optval)
      var.value[:] = optsol
      println("Mosek returns NaN!")
      println("gapopt = $gap")
      break
    end

    lastt = Problem(q.head, q.objective, q.constraints)

    inacct = t.value                             # the inaccurate value of t
    acct = maximum(evaluate.(constraints))       # the accurate value of t
    println("inaccuracy = $(inacct - acct)")
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
      println("SubOpt")
    end

    if lower < cur_low                           # update the lower bound
      low_sol = var.value[:]
      haha = copy(subgradients)
      qq = Problem(q.head, q.objective, q.constraints)
    else
      println("lowstab")
    end
    lower = max(cur_low, lower)
    println("upper = $upper")
    println("lower = $lower")
    gap = upper - lower
    if gap < TOL1
      var.value[:] = optsol
      println("gapopt = $gap")
      break
    else
      println("gapsub = $gap")
    end
    if upper < cur_upper
      rep += 1
      if length(permset) >= maxperm
        var.value[:] = optsol
        println("reached the maximum iteration at a single point")
        println("gapopt = $gap")
        break
      end
      if rep == checkpoint                         # check if the solver is stuck
        while pernum > maxperm
          println("reduce permnum!")
          lsTOL *= 0.1
          (c, pernum, perm) = greedy_rand(f.func, optsol, lsTOL)
        end
      end
      b = greedy(f.func, var.value[:])
      while in(b, subgradients)
        println("b included!")
        (b, pernum1, perm) = greedy_rand(f.func, var.value[:], lsTOL)
      end
      var.value[:] = optsol
      while length(permset) >= pernum
        println("length(permset) = $(length(permset))")
        println("pernum = $pernum")
        println("enlarge permnum!")
        lsTOL *= 1.1
        (a, pernum, perm) = greedy_rand(f.func, optsol, lsTOL)
      end
      (a, pernum, perm) = greedy_rand(f.func, optsol, lsTOL)
      while in(a, subgradients)
        println("a included!")
        if !in(perm, permset)
          push!(permset, perm)
        end
        (a, pernum) = greedy_rand(f.func, optsol, lsTOL)
      end
      push!(permset, perm)
      println(rep)
    else                                           # update the upper bound
      if max_rep < rep
        max_rep = copy(rep)
      end
      rep = 1                                      # reset repetition count
      lsTOL = 1e-3                                 # reset the level set tolerance
      (a, pernum, perm) = greedy_rand(f.func, optsol, lsTOL)
      a = greedy(f.func, optsol)
      permset = [sortperm(var.value[:], rev = true)]
      println("rep = $rep")
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
    println("cos_num = $cos_num")
    iters += 1
  end
  println("iters = $iters")
  println("max_rep = $max_rep")
  println("cos_num = $cos_num")
  return optsol, qq, lastt, low_sol
end

function cutting_plane(p::Problem, f::LovaszExtAbsAtom)

  lsTOL = 1e-3

  if length(get_v(p.objective + f)) != 1
    error("Cannot optimize functions with multiple variables.")
  end
  q = Problem(p.head, p.objective, p.constraints)
  # obj = p.objective + f
  conslength = length(p.constraints)

  # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = 1e-3))
  # solve!(q, SCSSolver(verbose=false))
  # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
  solve!(q, MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-6))

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
  haha = []

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

  for i = 1:maxiters

    println(length(permset))

    # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = TOL^2), warmstart = true)
    # solve!(q, SCSSolver(verbose=false, suppress_warnings = true), warmstart = true)
    # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
    solve!(q, MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-6))

    if isnan(q.optval)
      var.value[:] = optsol
      println("Mosek returns NaN!")
      println("gapopt = $gap")
      break
    end

    lastt = Problem(q.head, q.objective, q.constraints)

    inacct = t.value                             # the inaccurate value of t
    acct = maximum(evaluate.(constraints))       # the accurate value of t
    println("inaccuracy = $(inacct - acct)")
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
      println("SubOpt")
    end

    if lower < cur_low                           # update the lower bound
      low_sol = var.value[:]
      haha = copy(subgradients)
      qq = Problem(q.head, q.objective, q.constraints)
    else
      println("lowstab")
    end
    lower = max(cur_low, lower)
    println("upper = $upper")
    println("lower = $lower")
    gap = upper - lower
    if gap < TOL1
      var.value[:] = optsol
      println("gapopt = $gap")
      break
    else
      println("gapsub = $gap")
    end
    if upper < cur_upper
      if length(permset) >= maxperm
        var.value[:] = optsol
        println("reached the maximum iteration at a single point")
        println("gapopt = $gap")
        break
      end
      if rep == checkpoint                         # check if the solver is stuck
        while pernum > maxperm
          println("reduce permnum!")
          lsTOL *= 0.1
          (c, pernum, perm) = greedy_rand(f.func, abs.(optsol), lsTOL)
        end
      end
      rep += 1
      b = greedy(f.func, abs.(var.value[:]))
      while in(b, subgradients)
        println("b included!")
        (b, pernum1, perm) = greedy_rand(f.func, abs.(var.value[:]), lsTOL)
      end
      var.value[:] = optsol
      while length(permset) >= pernum
        println("length(permset) = $(length(permset))")
        println("pernum = $pernum")
        println("enlarge permnum!")
        lsTOL *= 1.1
        (a, pernum, perm) = greedy_rand(f.func, abs.(optsol), lsTOL)
      end
      (a, pernum, perm) = greedy_rand(f.func, abs.(optsol), lsTOL)
      while in(a, subgradients)
        println("a included!")
        if !in(perm, permset)
          push!(permset, perm)
        end
        (a, pernum) = greedy_rand(f.func, abs.(optsol), lsTOL)
      end
      push!(permset, perm)
      println(rep)
    else                                           # update the upper bound
      rep = 1                                      # reset repetition count
      lsTOL = 1e-3                                 # reset the level set tolerance
      (a, pernum, perm) = greedy_rand(f.func, abs.(optsol), lsTOL)
      a = greedy(f.func, abs.(optsol))
      permset = [sortperm(var.value[:], rev = true)]
      println(rep)
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
    println("cos_num = $cos_num")
    iters += 1
  end
  println("iters = $iters")
  println("cos_num = $cos_num")
  return optsol, qq, lastt, low_sol
end

function cutting_plane1(p::Problem, f::LovaszExtAtom)

  lsTOL = 1e-3

  if length(get_v(p.objective + f)) != 1
    error("Cannot optimize functions with multiple variables.")
  end
  q = Problem(p.head, p.objective, p.constraints)
  # obj = p.objective + f
  conslength = length(p.constraints)

  # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = 1e-3))
  # solve!(q, SCSSolver(verbose=false))
  solve!(q, ECOSSolver(verbose=false, abstol = 1e-4))
  # solve!(q, MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-6))

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
  haha = []

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

  for i = 1:maxiters

    println(length(permset))

    # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = TOL^2), warmstart = true)
    # solve!(q, SCSSolver(verbose=false, suppress_warnings = true), warmstart = true)
    solve!(q, ECOSSolver(verbose=false, abstol = 1e-4))
    # solve!(q, MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-6))

    if isnan(q.optval)
      var.value[:] = optsol
      println("Mosek returns NaN!")
      println("gapopt = $gap")
      break
    end

    lastt = Problem(q.head, q.objective, q.constraints)

    inacct = t.value                             # the inaccurate value of t
    acct = maximum(evaluate.(constraints))       # the accurate value of t
    println("inaccuracy = $(inacct - acct)")
    cur_sol = var.value[:]
    println("q.optval = $(q.optval)")
    cur_upper = q.optval - inacct + evaluate(f)[1]
    upper = min(upper, cur_upper)

    if q.status == :Optimal
      cur_low = q.optval - inacct + acct
      # cur_low = q.optval
    else
      cur_low = solvedualECOS!(q) - inacct + acct
      println("SubOpt")
    end

    if lower < cur_low                           # update the lower bound
      low_sol = var.value[:]
      haha = copy(subgradients)
      qq = Problem(q.head, q.objective, q.constraints)
    else
      println("lowstab")
    end
    lower = max(cur_low, lower)
    println("upper = $upper")
    println("lower = $lower")
    gap = upper - lower
    if gap < TOL1
      var.value[:] = optsol
      println("gapopt = $gap")
      break
    else
      println("gapsub = $gap")
    end
    if upper < cur_upper
      if length(permset) >= maxperm
        var.value[:] = optsol
        println("reached the maximum iteration at a single point")
        println("gapopt = $gap")
        break
      end
      if rep == checkpoint                         # check if the solver is stuck
        while pernum > maxperm
          println("reduce permnum!")
          lsTOL *= 0.1
          (c, pernum, perm) = greedy_rand(f.func, optsol, lsTOL)
        end
      end
      rep += 1
      b = greedy(f.func, var.value[:])
      while in(b, subgradients)
        println("b included!")
        (b, pernum1, perm) = greedy_rand(f.func, var.value[:], lsTOL)
      end
      var.value[:] = optsol
      while length(permset) >= pernum
        println("length(permset) = $(length(permset))")
        println("pernum = $pernum")
        println("enlarge permnum!")
        lsTOL *= 1.1
        (a, pernum, perm) = greedy_rand(f.func, optsol, lsTOL)
      end
      (a, pernum, perm) = greedy_rand(f.func, optsol, lsTOL)
      while in(a, subgradients)
        println("a included!")
        if !in(perm, permset)
          push!(permset, perm)
        end
        (a, pernum) = greedy_rand(f.func, optsol, lsTOL)
      end
      push!(permset, perm)
      println(rep)
    else                                           # update the upper bound
      rep = 1                                      # reset repetition count
      optsol = var.value[:]
      lsTOL = 1e-3                                 # reset the level set tolerance
      a = greedy(f.func, optsol)
      (a, pernum, perm) = greedy_rand(f.func, optsol, lsTOL)
      permset = [sortperm(var.value[:], rev = true)]
      println(rep)
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
    println("cos_num = $cos_num")
    iters += 1
  end
  println("iters = $iters")
  println("cos_num = $cos_num")
  return optsol, qq, lastt, low_sol
end

# function cutting_plane1(p::Problem, f::LovaszExtAtom)
#
#   lsTOL = 1e-3
#
#   if length(get_v(p.objective)) > 1
#     error("Cannot optimize functions with multiple variables.")
#   end
#   q = Problem(p.head, p.objective, p.constraints)
#   obj = p.objective + f
#   conslength = length(p.constraints)
#
#   # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = 1e-6))
#   # solve!(q, SCSSolver(verbose=false))
#   solve!(q, ECOSSolver(verbose=false, abstol = 1e-6, feastol = 1e-6))
#   # solve!(q, MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-6))
#
#   # Reset the primal variable for warmstart
#   push!(q.solution.primal, q.solution.primal[end])
#   q.solution.primal[end - 1] = 0
#
#   t = Variable(1)
#   q.objective += t
#   var = get_v(q.objective)[1]
#
#   n = size(var)[1]
#
#   a = greedy(f.func, var.value[:])
#   push!(q.constraints, dot(a, var) <= t)
#   constraints = []
#   push!(constraints, q.constraints[end].lhs)       # the array of the affine functions
#   subgradients = [a]
#   haha = []
#
#   cur_upper = 0.0
#   cur_low = 0.0
#   upper = Inf
#   lower = -Inf
#   optsol = var.value[:]
#   pernum = 1
#   permset = [sortperm(optsol, rev = true)]         # orderings explored
#   low_sol = zeros(length(var.value[:]))
#   iters = 1
#   rep = 0                                          # the repetition count of a trial point
#   cos_num = 0.0
#   gap = 0.0
#
#   lastt = Problem(q.head, q.objective, q.constraints)
#
#   qq = Problem(q.head, q.objective, q.constraints)
#
#   (objvec, A, dualvec, cones, var_to_ranges, vartypes, conic_constraints) = conic_problem(q)
#
#   for i = 1:maxiters
#
#     println(length(permset))
#
#     # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = 1e-6), warmstart = true)
#     # solve!(q, SCSSolver(verbose=false, suppress_warnings = true), warmstart = true)
#     solve!(q, ECOSSolver(verbose=false, abstol = 1e-6, feastol = 1e-6))
#     # solve!(q, MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-6))
#
#     if isnan(q.optval)
#       var.value[:] = optsol
#       println("Mosek returns NaN!")
#       println("gapopt = $gap")
#       break
#     end
#
#     lastt = Problem(q.head, q.objective, q.constraints)
#
#     if q.status == :Optimal
#       inacct = t.value                            # the inaccurate value of t
#       acct = maximum(evaluate.(constraints))      # the accurate value of t
#       cur_sol = var.value[:]
#       cur_upper = q.optval - inacct + evaluate(f)[1]
#       upper = min(upper, cur_upper)
#       cur_low = q.optval - inacct + acct
#
#       if lower < cur_low                           # update the lower bound
#         low_sol = var.value[:]
#         haha = copy(subgradients)
#         qq = Problem(q.head, q.objective, q.constraints)
#       else
#         println("lowstab")
#       end
#       lower = max(cur_low, lower)
#     else
#       println("SubOpt")
#     end
#     println("upper = $upper")
#     println("lower = $lower")
#     gap = upper - lower
#     if gap < TOL1
#       var.value[:] = optsol
#       println("gapopt = $gap")
#       break
#     else
#       println("gapsub = $gap")
#     end
#     if upper < cur_upper
#       if length(permset) >= maxperm
#         var.value[:] = optsol
#         println("reached the maximum iteration at a single point")
#         println("gapopt = $gap")
#         break
#       end
#       if rep == checkpoint                         # check if the solver is stuck
#         while pernum > maxperm
#           println("reduce permnum!")
#           lsTOL *= 0.1
#           (c, pernum, perm) = greedy_rand(f.func, optsol, lsTOL)
#         end
#       end
#       rep += 1
#       b = greedy(f.func, var.value[:])
#       while in(b, subgradients)
#         println("b included!")
#         (b, pernum1, perm) = greedy_rand(f.func, var.value[:], lsTOL)
#       end
#       var.value[:] = optsol
#       while length(permset) >= pernum
#         println("length(permset = $(length(permset))")
#         println("pernum = $pernum")
#         println("enlarge permnum!")
#         lsTOL *= 1.1
#         (a, pernum, perm) = greedy_rand(f.func, optsol, lsTOL)
#       end
#       (a, pernum, perm) = greedy_rand(f.func, optsol, lsTOL)
#       while in(a, subgradients)
#         println("a included!")
#         if !in(perm, permset)
#           push!(permset, perm)
#         end
#         (a, pernum) = greedy_rand(f.func, optsol, lsTOL)
#       end
#       push!(permset, perm)
#       println(rep)
#     else                                           # update the upper bound
#       rep = 1                                      # reset repetition count
#       optsol = var.value[:]
#       lsTOL = 1e-3                                 # reset the level set tolerance
#       a = greedy(f.func, optsol)
#       (a, pernum, perm) = greedy_rand(f.func, optsol, lsTOL)
#       permset = [sortperm(var.value[:], rev = true)]
#       println(rep)
#     end
#
#     cos_num = length(constraints)
#     if cos_num > n
#       prod = zeros(cos_num)
#       for i = 1:cos_num
#         prod[i] = dot(optsol, subgradients[i])
#       end
#       if cos_num == (n + 2)
#         minn = indmin(prod)
#         minn = indmin(prod)
#         q.constraints[conslength + minn] = q.constraints[end]
#         constraints[minn] = constraints[end]
#         subgradients[minn] = subgradients[end]
#         prod = prod[1:end-1]
#         q.constraints = q.constraints[1:end-1]
#         constraints = constraints[1:end-1]
#         subgradients = subgradients[1:end-1]
#       end
#       if rep == 1
#         minn = indmin(prod)
#         q.constraints[conslength + minn] = (dot(a, var) <= t)
#         constraints[minn] = q.constraints[conslength + minn].lhs
#         subgradients[minn] = a
#       else                                       # b is added
#         minn = indmin(prod)
#         q.constraints[conslength + minn] = (dot(a, var) <= t)
#         constraints[minn] = q.constraints[conslength + minn].lhs
#         subgradients[minn] = a
#         maxx = indmin(prod)
#         prod[minn] = prod[maxx]
#         minn1 = indmin(prod)
#         q.constraints[conslength + minn1] = (dot(b, var) <= t)
#         constraints[minn1] = q.constraints[conslength + minn1].lhs
#         subgradients[minn1] = b
#       end
#     else
#       if rep > 1
#         push!(q.constraints, dot(b, var) <= t)
#         push!(constraints, q.constraints[end].lhs)
#         push!(subgradients, b)
#       end
#       push!(q.constraints, dot(a, var) <= t)
#       push!(constraints, q.constraints[end].lhs)
#       push!(subgradients, a)
#     end
#     println("cos_num = $cos_num")
#     iters += 1
#   end
#   println("iters = $iters")
#   println("cos_num = $cos_num")
#   return optsol, qq, lastt, low_sol
# end
#
# function cutting_planeMosek(p::Problem, f::LovaszExtAtom)
#
#   lsTOL = 1e-3
#
#   if length(get_v(p.objective)) > 1
#     error("Cannot optimize functions with multiple variables.")
#   end
#   q = Problem(p.head, p.objective, p.constraints)
#   obj = p.objective + f
#   conslength = length(p.constraints)
#
#   solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = 1e-3))
#   # solve!(q, SCSSolver(verbose=false))
#   # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
#   # solve!(q, MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-6))
#
#   # Reset the primal variable for warmstart
#   push!(q.solution.primal, q.solution.primal[end])
#   q.solution.primal[end - 1] = 0
#
#   t = Variable(1)
#   q.objective += t
#   var = get_v(q.objective.children[1])[1]
#   n = size(var)[1]
#   a = greedy(f.func, var.value[:])
#   push!(q.constraints, dot(a, var) <= t)
#   constraints = []
#   push!(constraints, q.constraints[end].lhs)       # the array of the affine functions
#   subgradients = [a]
#   haha = []
#
#   cur_upper = 0.0
#   cur_low = 0.0
#   upper = Inf
#   lower = -Inf
#   optsol = var.value[:]
#   pernum = 1
#   permset = [sortperm(optsol, rev = true)]         # orderings explored
#   low_sol = zeros(length(var.value[:]))
#   iters = 1
#   rep = 0                                          # the repetition count of a trial point
#   cos_num = 0.0
#   gap = 0.0
#
#   lastt = Problem(q.head, q.objective, q.constraints)
#
#   qq = Problem(q.head, q.objective, q.constraints)
#
#   for i = 1:maxiters
#
#     println(length(permset))
#
#     solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = TOL^2), warmstart = true)
#     # solve!(q, SCSSolver(verbose=false, suppress_warnings = true), warmstart = true)
#     # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
#     # solve!(q, MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-6))
#
#     if isnan(q.optval)
#       var.value[:] = optsol
#       println("Mosek returns NaN!")
#       println("gapopt = $gap")
#       break
#     end
#
#     lastt = Problem(q.head, q.objective, q.constraints)
#
#     if q.status == :Optimal
#       cur_sol = var.value[:]
#       cur_upper = q.optval - t.value + evaluate(f)[1]
#       upper = min(upper, cur_upper)
#       cur_low = q.optval
#       # cur_low = q.optval
#       if lower < cur_low                           # update the lower bound
#         low_sol = var.value[:]
#         haha = copy(subgradients)
#         qq = Problem(q.head, q.objective, q.constraints)
#       else
#         println("lowstab")
#       end
#       lower = max(cur_low, lower)
#     else
#       println("SubOpt")
#     end
#     println("upper = $upper")
#     println("lower = $lower")
#     gap = upper - lower
#     if gap < TOL1
#       var.value[:] = optsol
#       println("gapopt = $gap")
#       break
#     else
#       println("gapsub = $gap")
#     end
#     if upper < cur_upper
#       if length(permset) >= maxperm
#         var.value[:] = optsol
#         println("reached the maximum iteration at a single point")
#         println("gapopt = $gap")
#         break
#       end
#       if rep == checkpoint                         # check if the solver is stuck
#         while pernum > maxperm
#           println("reduce permnum!")
#           lsTOL *= 0.1
#           (c, pernum, perm) = greedy_rand(f.func, optsol, lsTOL)
#         end
#       end
#       rep += 1
#       b = greedy(f.func, var.value[:])
#       while in(b, subgradients)
#         println("b included!")
#         (b, pernum1, perm) = greedy_rand(f.func, var.value[:], lsTOL)
#       end
#       var.value[:] = optsol
#       while length(permset) >= pernum
#         println("length(permset = $(length(permset))")
#         println("pernum = $pernum")
#         println("enlarge permnum!")
#         lsTOL *= 1.1
#         (a, pernum, perm) = greedy_rand(f.func, optsol, lsTOL)
#       end
#       (a, pernum, perm) = greedy_rand(f.func, optsol, lsTOL)
#       while in(a, subgradients)
#         println("a included!")
#         if !in(perm, permset)
#           push!(permset, perm)
#         end
#         (a, pernum) = greedy_rand(f.func, optsol, lsTOL)
#       end
#       push!(permset, perm)
#       println(rep)
#     else                                           # update the upper bound
#       rep = 1                                      # reset repetition count
#       optsol = var.value[:]
#       lsTOL = 1e-3                                 # reset the level set tolerance
#       a = greedy(f.func, optsol)
#       (a, pernum, perm) = greedy_rand(f.func, optsol, lsTOL)
#       permset = [sortperm(var.value[:], rev = true)]
#       println(rep)
#     end
#
#     cos_num = length(constraints)
#     if cos_num > n
#       prod = zeros(cos_num)
#       for i = 1:cos_num
#         prod[i] = dot(optsol, subgradients[i])
#       end
#       if cos_num == (n + 2)
#         minn = indmin(prod)
#         minn = indmin(prod)
#         q.constraints[conslength + minn] = q.constraints[end]
#         constraints[minn] = constraints[end]
#         subgradients[minn] = subgradients[end]
#         prod = prod[1:end-1]
#         q.constraints = q.constraints[1:end-1]
#         constraints = constraints[1:end-1]
#         subgradients = subgradients[1:end-1]
#       end
#       if rep == 1
#         minn = indmin(prod)
#         q.constraints[conslength + minn] = (dot(a, var) <= t)
#         constraints[minn] = q.constraints[conslength + minn].lhs
#         subgradients[minn] = a
#       else                                       # b is added
#         minn = indmin(prod)
#         q.constraints[conslength + minn] = (dot(a, var) <= t)
#         constraints[minn] = q.constraints[conslength + minn].lhs
#         subgradients[minn] = a
#         maxx = indmin(prod)
#         prod[minn] = prod[maxx]
#         minn1 = indmin(prod)
#         q.constraints[conslength + minn1] = (dot(b, var) <= t)
#         constraints[minn1] = q.constraints[conslength + minn1].lhs
#         subgradients[minn1] = b
#       end
#     else
#       if rep > 1
#         push!(q.constraints, dot(b, var) <= t)
#         push!(constraints, q.constraints[end].lhs)
#         push!(subgradients, b)
#       end
#       push!(q.constraints, dot(a, var) <= t)
#       push!(constraints, q.constraints[end].lhs)
#       push!(subgradients, a)
#     end
#     println("cos_num = $cos_num")
#     iters += 1
#   end
#   println("iters = $iters")
#   println("cos_num = $cos_num")
#   return optsol, qq, lastt, low_sol
# end

# function cutting_plane(p::Problem, f::LovaszExtAtom)
#   if length(get_v(p.objective)) > 1
#     error("Cannot optimize functions with multiple variables.")
#   end
#   q = Problem(p.head, p.objective, p.constraints)
#   obj = p.objective + f
#
#   # solve!(q, SCSSolver(verbose=false))
#   # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
#   # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = 1e-3))
#   solve!(q, MosekSolver(MSK_IPAR_LOG = 0))
#
#   # Reset the primal variable for warmstart
#   push!(q.solution.primal, q.solution.primal[end])
#   q.solution.primal[end - 1] = 0
#
#   t = Variable(1)
#   q.objective += t
#   var = get_v(q.objective.children[1])[1]
#   a = greedy_rand(f.func, var.value[:])
#   push!(q.constraints, dot(a, var) <= t)
#   constraints = []
#   push!(constraints, q.constraints[end].lhs)       # The array of the affine functions
#   subgradients = [a]
#
#   cur_upper = 0.0
#   cur_low = 0.0
#   upper = evaluate(obj)[1]
#   lower = 0.0
#   optsol = var.value[:]
#   low_sol = var.value[:]
#   iters = 0
#   rep = 0                                          # The repetition count of a trial point
#
#   qq = Problem(q.head, q.objective, q.constraints)
#
#   for i = 1:maxiters
#     # println()
#     # solve!(q, SCSSolver(verbose=false, suppress_warnings = true), warmstart = true)
#
#     # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
#
#     # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = TOL^2), warmstart = true)
#
#     solve!(q, MosekSolver(MSK_IPAR_LOG = 0))
#
#     cur_sol = copy(var.value[:])
#
#     # solve!(q, MosekSolver(MSK_IPAR_LOG = 0))
#     # t.value = maximum(evaluate.(constraints))
#     # lowermos = evaluate(q.objective)[1]
#     # uppermos = evaluate(obj)[1]
#     #
#     # var.value[:] = copy(cur_sol)
#
#     # cur_sol = copy(var.value[:])
#
#     cur_upper = evaluate(obj)[1]
#     if upper > cur_upper
#       upper = cur_upper
#       optsol = var.value[:]
#     end
#     t.value = maximum(evaluate.(constraints))      # Adjust the value of t
#     if q.status == :Optimal
#       cur_low = evaluate(q.objective)[1]
#       println("lower_improv = $(cur_low - lower)")
#       if lower < cur_low
#         low_sol = var.value[:]
#       end
#       lower = max(cur_low, lower)
#     else
#       println("SubOpt")
#     end
#
#     # println("mosup - curup = $(uppermos - cur_upper)")
#     # println("curlow - moslowe = $(cur_low - lowermos)")
#
#     println("upper = $upper")
#     println("lower = $lower")
#     gap = upper - lower
#     if gap < TOL
#       println("gap = $gap")
#       break
#     else
#       println("gap = $gap")
#     end
#     if upper < cur_upper
#       rep += 1
#       b = greedy_rand(f.func, var.value[:])
#       while in(b, subgradients)
#         b = greedy_rand(f.func, var.value[:])
#       end
#       prod = zeros(length(subgradients))
#       for i = 1:length(subgradients)
#         prod[i] = dot(var.value[:], subgradients[i])
#       end
#       push!(q.constraints, dot(b, var) <= t)
#       push!(constraints, q.constraints[end].lhs)
#       push!(subgradients, b)
#
#       # if maxx - minn > 1e-2
#       #   q.constraints[minn] = (dot(b, var) <= t)
#       #   constraints[minn] = q.constraints[minn].lhs
#       #   subgradients[minn] = b
#       # else
#       #   push!(q.constraints, dot(b, var) <= t)
#       #   push!(constraints, q.constraints[end].lhs)
#       #   push!(subgradients, b)
#       # end
#
#       var.value = copy(optsol)
#       if rep == 1
#         a = greedy(f.func, var.value[:])
#       else
#         a = greedy_rand(f.func, var.value[:])
#         while in(a, subgradients)
#           println("included~")
#           a = greedy_rand(f.func, var.value[:])
#         end
#
#         prod = zeros(length(subgradients))
#         for i = 1:length(subgradients)
#           prod[i] = dot(var.value[:], subgradients[i])
#         end
#         minn = indmin(prod)
#         maxx = indmax(prod)
#
#         if maxx - minn > 1e-2
#           q.constraints[minn] = (dot(a, var) <= t)
#           constraints[minn] = q.constraints[minn].lhs
#           subgradients[minn] = a
#         else
#           push!(q.constraints, dot(a, var) <= t)
#           push!(constraints, q.constraints[end].lhs)
#           push!(subgradients, a)
#         end
#
#       end
#       println(rep)
#     else
#       rep = 0                                      # Reset repetition count
#       optsol = var.value[:]
#       a = greedy_rand(f.func, var.value[:])
#       while in(a, subgradients)
#         a = greedy_rand(f.func, var.value[:])
#       end
#       println(rep)
#     end
#
#     prod = zeros(length(subgradients))
#     for i = 1:length(subgradients)
#       prod[i] = dot(var.value[:], subgradients[i])
#     end
#     minn = indmin(prod)
#     maxx = indmax(prod)
#
#     if maxx - minn > 1e-2
#       q.constraints[minn] = (dot(a, var) <= t)
#       constraints[minn] = q.constraints[minn].lhs
#       subgradients[minn] = a
#     else
#       push!(q.constraints, dot(a, var) <= t)
#       push!(constraints, q.constraints[end].lhs)
#       push!(subgradients, a)
#     end
#
#     iters += 1
#   end
#   println(iters)
#   println("length = $(length(subgradients))")
#   upper_low = zeros(length(optsol), 2)
#   upper_low[:, 1] = optsol
#   upper_low[:, 2] = low_sol
#   return upper_low
# end

# function cutting_plane(p::Problem, f::LovaszExtAtom)
#   if length(get_v(p.objective)) > 1
#     error("Cannot optimize functions with multiple variables.")
#   end
#   q = Problem(p.head, p.objective, p.constraints)
#   obj = p.objective + f
#
#   # solve!(q, SCSSolver(verbose=false))
#   # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
#   # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = 1e-3))
#   solve!(q, MosekSolver(MSK_IPAR_LOG = 0))
#
#   # Reset the primal variable for warmstart
#   push!(q.solution.primal, q.solution.primal[end])
#   q.solution.primal[end - 1] = 0
#
#   t = Variable(1)
#   q.objective += t
#   var = get_v(q.objective.children[1])[1]
#   a = greedy_rand(f.func, var.value[:])
#   push!(q.constraints, dot(a, var) <= t)
#   constraints = []
#   push!(constraints, q.constraints[end].lhs)       # The array of the affine functions
#   subgradients = [a]
#
#   cur_upper = 0.0
#   cur_low = 0.0
#   upper = evaluate(obj)[1]
#   lower = 0.0
#   optsol = var.value[:]
#   low_sol = var.value[:]
#   iters = 0
#   rep = 0                                          # The repetition count of a trial point
#
#   qq = Problem(q.head, q.objective, q.constraints)
#
#   for i = 1:maxiters
#     # println()
#     # solve!(q, SCSSolver(verbose=false, suppress_warnings = true), warmstart = true)
#
#     # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
#
#     # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = TOL^2), warmstart = true)
#
#     rep = 0
#
#     solve!(q, MosekSolver(MSK_IPAR_LOG = 0))
#
#     if q.status == :Optimal
#       cur_low = evaluate(q.objective)[1]
#       # println("lower_improv = $(cur_low - lower)")
#       # println("low_diff = ($(cur_low - evaluate(obj)[1]))")
#       if lower < cur_low
#         low_sol = var.value[:]
#       end
#       lower = max(cur_low, lower)
#     else
#       println("SubOpt")
#     end
#
#     cur_upper = evaluate(obj)[1]
#     upper = min(upper, cur_upper)
#     gap = upper - lower
#
#     if gap < TOL
#       println("upper = $upper")
#       println("lower = $lower")
#       println("gap = $gap")
#       break
#     end
#
#     cur_subgradients = [a]
#
#     while upper <= cur_upper
#       rep += 1
#       gap = upper - lower
#       if gap < TOL
#         println("upper = $upper")
#         println("lower = $lower")
#         println("gap = $gap")
#         break
#       end
#       while in(a, cur_subgradients)
#         a = greedy_rand(f.func, optsol)
#       end
#       q.constraints[end] = (dot(a, var) <= t)
#       constraints[end] = q.constraints[end].lhs
#       subgradients[end] = a
#       push!(cur_subgradients, a)
#       solve!(q, MosekSolver(MSK_IPAR_LOG = 0))
#       if q.status == :Optimal
#         t.value = maximum(evaluate.(constraints))      # Adjust the value of t
#         cur_low = evaluate(q.objective)[1]
#         # println("lower_improv = $(cur_low - lower)")
#         # println("low_diff = ($(cur_low - evaluate(obj)[1]))")
#         if lower < cur_low
#           low_sol = var.value[:]
#         end
#         lower = max(cur_low, lower)
#       else
#         println("SubOpt")
#       end
#       cur_upper = evaluate(obj)[1]
#       println("upper = $upper")
#       println("lower = $lower")
#       println("gap = $gap")
#       println("rep")
#     end
#
#     cur_sol = copy(var.value[:])
#
#     upper = min(cur_upper, upper)
#
#     # solve!(q, MosekSolver(MSK_IPAR_LOG = 0))
#     # t.value = maximum(evaluate.constraints)
#     # lowermos = evaluate(q.objective)[1]
#     # uppermos = evaluate(obj)[1]
#     #
#     # var.value[:] = copy(cur_sol)
#
#     # cur_sol = copy(var.value[:])
#
#     # if q.status == :Optimal
#     #   t.value = maximum(evaluate.(constraints))      # Adjust the value of t
#     #   cur_low = evaluate(q.objective)[1]
#     #   # println("lower_improv = $(cur_low - lower)")
#     #   # println("low_diff = ($(cur_low - evaluate(obj)[1]))")
#     #   if lower < cur_low
#     #     low_sol = var.value[:]
#     #   end
#     #   lower = max(cur_low, lower)
#     # else
#     #   println("SubOpt")
#     # end
#
#     # println("Mosup - up = $(upper - uppermos)")
#     # println("lower - Moslower = $(lower - lowermos)")
#
#     # println("upper = $upper")
#     # println("lower = $lower")
#     # gap = upper - lower
#     # if gap < TOL
#     #   println("gap = $gap")
#     #   break
#     # else
#     #   println("gap = $gap")
#     # end
#     # if upper < cur_upper
#     #   rep += 1
#     #   b = greedy_rand(f.func, var.value[:])
#     #   while in(b, subgradients)
#     #     b = greedy_rand(f.func, var.value[:])
#     #   end
#     #   push!(q.constraints, dot(b, var) <= t)
#     #   push!(constraints, q.constraints[end].lhs)
#     #   push!(subgradients, b)
#     #   var.value = copy(optsol)
#     #   if rep == 1
#     #     a = greedy(f.func, var.value[:])
#     #   else
#     #     a = greedy_rand(f.func, var.value[:])
#     #     while in(a, subgradients)
#     #       a = greedy_rand(f.func, var.value[:])
#     #     end
#     #   end
#     #   println(rep)
#     # else
#     #   rep = 0                                      # Reset repetition count
#     #   optsol = var.value[:]
#     #   a = greedy_rand(f.func, var.value[:])
#     #   while in(a, subgradients)
#     #     a = greedy_rand(f.func, var.value[:])
#     #   end
#     #   println(rep)
#     # end
#     optsol = var.value[:]
#     a = greedy_rand(f.func, optsol)
#     push!(q.constraints, dot(a, var) <= t)
#     push!(constraints, q.constraints[end].lhs)
#     push!(subgradients, a)
#     iters += 1
#   end
#   println(iters)
#   upper_low = zeros(length(optsol), 2)
#   upper_low[:, 1] = optsol
#   upper_low[:, 2] = low_sol
#   return upper_low
# end

# function cutting_planeECOS(p::Problem, f::LovaszExtAtom)
#   if length(get_v(p.objective)) > 1
#     error("Cannot optimize functions with multiple variables.")
#   end
#   q = Problem(p.head, p.objective, p.constraints)
#   obj = p.objective + f
#
#   # solve!(q, SCSSolver(verbose=false))
#   # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
#   solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = 1e-3))
#   # solve!(q, MosekSolver(MSK_IPAR_LOG = 0))
#
#   # Reset the primal variable for warmstart
#   push!(q.solution.primal, q.solution.primal[end])
#   q.solution.primal[end - 1] = 0
#
#   t = Variable(1)
#   q.objective += t
#   var = get_v(q.objective.children[1])[1]
#   a = greedy_rand(f.func, var.value[:])
#   push!(q.constraints, dot(a, var) <= t)
#   constraints = []
#   push!(constraints, q.constraints[end].lhs)       # The array of the affine functions
#   subgradients = [a]
#
#   cur_upper = 0.0
#   cur_low = 0.0
#   upper = evaluate(obj)[1]
#   lower = 0.0
#   optsol = var.value[:]
#   low_sol = var.value[:]
#   iters = 0
#   rep = 0                                          # The repetition count of a trial point
#
#   for i = 1:maxiters
#     # println()
#     # solve!(q, SCSSolver(verbose=false, suppress_warnings = true), warmstart = true)
#
#     # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
#
#     solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = TOL^2), warmstart = true)
#
#     # solve!(q, MosekSolver(MSK_IPAR_LOG = 0))
#
#     cur_sol = copy(var.value[:])
#
#     # solve!(q, MosekSolver(MSK_IPAR_LOG = 0))
#     # t.value = maximum(evaluate.constraints)
#     # lowermos = evaluate(q.objective)[1]
#     # uppermos = evaluate(obj)[1]
#     #
#     # var.value[:] = copy(cur_sol)
#
#     # cur_sol = copy(var.value[:])
#
#     cur_upper = evaluate(obj)[1]
#     if upper > cur_upper
#       upper = cur_upper
#       optsol = var.value[:]
#     end
#     t.value = maximum(evaluate.(constraints))      # Adjust the value of t
#     if q.status == :Optimal
#       cur_low = evaluate(q.objective)[1]
#       # println("lower_improv = $(cur_low - lower)")
#       # println("low_diff = ($(cur_low - evaluate(obj)[1]))")
#       if lower < cur_low
#         low_sol = var.value[:]
#       end
#       lower = max(cur_low, lower)
#     else
#       println("SubOpt")
#     end
#
#     # println("Mosup - up = $(upper - uppermos)")
#     # println("lower - Moslower = $(lower - lowermos)")
#
#     println("upper = $upper")
#     println("lower = $lower")
#     gap = upper - lower
#     if gap < TOL
#       println("gap = $gap")
#       break
#     else
#       println("gap = $gap")
#     end
#     if upper < cur_upper
#       rep += 1
#       b = greedy_rand(f.func, var.value[:])
#       while in(b, subgradients)
#         b = greedy_rand(f.func, var.value[:])
#       end
#       push!(q.constraints, dot(b, var) <= t)
#       push!(constraints, q.constraints[end].lhs)
#       push!(subgradients, b)
#       var.value = copy(optsol)
#       if rep == 1
#         a = greedy(f.func, var.value[:])
#       else
#         a = greedy_rand(f.func, var.value[:])
#         while in(a, subgradients)
#           a = greedy_rand(f.func, var.value[:])
#         end
#       end
#       println(rep)
#     else
#       rep = 0                                      # Reset repetition count
#       optsol = var.value[:]
#       a = greedy_rand(f.func, optsol)
#       while in(a, subgradients)
#         a = greedy_rand(f.func, optsol)
#       end
#       println(rep)
#     end
#     push!(q.constraints, dot(a, var) <= t)
#     push!(constraints, q.constraints[end].lhs)
#     push!(subgradients, a)
#     iters += 1
#   end
#   println(iters)
#   upper_low = zeros(length(optsol), 2)
#   upper_low[:, 1] = optsol
#   upper_low[:, 2] = low_sol
#   return upper_low
# end

# function cutting_planeECOS(p::Problem, f::LovaszExtAtom)
#   if length(get_v(p.objective)) > 1
#     error("Cannot optimize functions with multiple variables.")
#   end
#   q = Problem(p.head, p.objective, p.constraints)
#   obj = p.objective + f
#
#   # solve!(q, SCSSolver(verbose=false))
#   solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
#   # solve!(q, GurobiSolver(OutputFlag=0), warmstart = true)
#
#   # Reset the primal variable for warmstart
#   push!(q.solution.primal, q.solution.primal[end])
#   q.solution.primal[end - 1] = 0
#
#   t = Variable(1)
#   q.objective += t
#   var = get_v(q.objective.children[1])[1]
#   a = greedy_rand(f.func, var.value[:])
#   push!(q.constraints, dot(a, var) <= t)
#   constraints = []
#   subgradients = []
#   push!(constraints, q.constraints[end].lhs)       # The array of the affine functions
#   push!(subgradients, a)
#
#   cur_upper = 0.0
#   cur_low = 0.0
#   upper = evaluate(obj)[1]
#   lower = 0.0
#   optsol = var.value[:]
#   here = var.value[:]
#   iters = 0
#   rep = 0                                          # The repetition count of a trial point
#
#   qq = Problem(q.head, q.objective, q.constraints)
#
#   for i = 1:maxiters
#     println()
#     # solve!(q, SCSSolver(verbose=false, suppress_warnings = true), warmstart = true)
#
#     solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
#     qq = Problem(q.head, q.objective, q.constraints)
#
#     # solve!(q, GurobiSolver(OutputFlag=0), warmstart = true)
#
#     here = var.value[:]
#
#     cur_upper = evaluate(obj)[1]
#     upper = min(upper, cur_upper)
#     t.value = maximum(evaluate.(constraints))      # Adjust the value of t
#     if q.status == :Optimal
#       cur_low = evaluate(q.objective)[1]
#       println("lower_improv = $(cur_low - lower)")
#       lower = max(cur_low, lower)
#     else
#       println("SubOpt")
#     end
#
#     println("upper = $upper")
#     println("lower = $lower")
#     gap = upper - lower
#     if gap < TOL
#       println("gap = $gap")
#       println(iters)
#       return var.value[:]
#       break
#     else
#       println("gap = $gap")
#     end
#     if upper < cur_upper
#       rep += 1
#       b = greedy_rand(f.func, var.value[:])
#       if !in(b, subgradients)
#         push!(q.constraints, dot(b, var) <= t)
#         push!(constraints, q.constraints[end].lhs)
#         push!(subgradients, b)
#       end
#       var.value = copy(optsol)
#       if rep == 1
#         a = greedy(f.func, var.value[:])
#       else
#         a = greedy_rand(f.func, var.value[:])
#       end
#       println(rep)
#     else
#       rep = 0                                      # Reset repetition count
#       upper = copy(cur_upper)
#       optsol = var.value[:]
#       a = greedy_rand(f.func, var.value[:])
#       println(rep)
#     end
#     if !in(a, subgradients)
#       push!(q.constraints, dot(a, var) <= t)
#       push!(constraints, q.constraints[end].lhs)
#       push!(subgradients, a)
#     end
#     iters += 1
#   end
#   println(iters)
#   return optsol, here
# end
#
# function cutting_plane1(p::Problem, f::LovaszExtAtom)
#   if length(get_v(p.objective)) > 1
#     error("Cannot optimize functions with multiple variables.")
#   end
#
#   # equa(x) = x > -1e-2
#
#   q = Problem(p.head, p.objective, p.constraints)
#   obj = p.objective + f
#   solve!(q, SCSSolver(verbose=false))
#   # solve!(q, ECOSSolver(verbose=false, absstol = 1e-6))
#   # solve!(q, GurobiSolver(OutputFlag=0))
#   # Reset the primal variable for warmstart
#   push!(q.solution.primal, q.solution.primal[end])
#   q.solution.primal[end - 1] = 0
#   q.solution.status = "not yet solved"
#   t = Variable(1)
#   q.objective += t
#   var = get_v(q.objective.children[1])[1]
#   a = greedy_rand(f.func, var.value[:])
#   push!(q.constraints, dot(a, var) <= t)
#   constraints = []
#   push!(constraints, q.constraints[end].lhs)       # The array of the affine functions
#   # upper = evaluate(obj)[1]
#   # optsol = var.value[:]
#   iters = 0
#   for i = 1:100
#     solve!(q, SCSSolver(verbose=false, suppress_warnings = true), warmstart = true)
#     # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
#     # solve!(q, GurobiSolver(OutputFlag=0), warmstart = true)
#     t.value = maximum(evaluate.(constraints))      # Adjusting the value of t
#     println("upper = $(evaluate(obj)[1])")
#     println("lower = $((evaluate(p.objective) + t.value)[1])")
#     # if upper > evaluate(obj)[1]
#     #   upper = evaluate(obj)[1]
#     #   optsol = var.value[:]
#     # end
#     # gap = upper - (evaluate(p.objective) + t.value)[1]
#     gap = evaluate(f, var.value[:, 1])[1] - t.value
#     if gap < TOL
#       println("gap = $gap")
#       println("diff = $(evaluate(f)[1] - t.value)")
#       return var.value[:]
#       break
#     else
#       println("gap = $gap")
#       println()
#       a = greedy_rand(f.func, var.value[:])
#       # equ_ind = find(equa, evaluate.(constraints) - t.value)
#       # println(evaluate.(constraints) - t.value)
#       # println(equ_ind)
#       # println()
#       # q.constraints = q.constraints[equ_ind]
#       # constraints = constraints[equ_ind]
#       push!(q.constraints, dot(a, var) <= t)
#       push!(constraints, q.constraints[end].lhs)
#     end
#     iters += 1
#   end
#   println(iters)
#   return var.value[:, 1]
# end
