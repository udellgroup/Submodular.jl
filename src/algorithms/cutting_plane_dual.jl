#############################################################################
# cutting_plane_dual.jl
# convert a problem of the DualModel to the PrimalModel and solve it with the
# cutting plane method.
# Conic problems having the form
# minimize c'*x
# st       b - Ax ∈ K, where K is a cone,
# their dual problems should be
# maximize -b'*y
# st       c + Aᵀy ∈ K^*, where K^* is the dual cone of K.
#############################################################################

export cutting_plane_dual

function cutting_plane_dual(prob::Problem,
                            p::BasePoly;
                            s::AbstractMathProgSolver = MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = abs_tol*1e-3),
                            abs_tol::Float64 = 1e-3,
                            max_iters::Int = 100,
                            max_rep::Int = 100,
                            lev_tol_def::Float64 = 1e-3,
                            lev_tol_up::Float64 = 1e-1,
                            lev_tol_low::Float64 = 1e-7,
                            shrink_point::Int = 20,
                            λ::Float64 = 1.5,
                            μ::Float64 = 0.7)

  lev_tol = lev_tol_def

  F = p.F
  if length(get_v(p.objective)) != 1
    error("Cannot optimize functions with multiple variables.")
  end
  # convert the problem to primal form
  c, A, b, dual_constr_cones, dual_var_cones = dual_to_primal(p)
  Ã = A'                                          # transpose of A, stored to implement faster row addition
  m = ConicModel(s)

  var = get_v(prob.objective)[1]
  n = size(var)[1]

  a = greedy(F, zeros(n))                         # first approximation of the Lovasz extension
  subgradients = [a]
  slope = zeros(size(A))
  push!(slope[1], -1)
  Ã =

  # initialization
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

  loadproblem!(m, vec(full(b)), -A', vec(full(c)), dual_constr_cones, dual_var_cones)
  optimize!(m)

  solution = try
    MathProgBase.getsolution(m)
  catch
    fill(NaN, numvar(m))
  end

  return dot(c, -solution)

  push!(q.constraints, dot(a, var) <= t)
  constraints = []
  push!(constraints, q.constraints[end].lhs)       # the array of the affine functions



  for i = 1:max_iters

    # println(length(permset))

    # solve!(q, GurobiSolver(OutputFlag=0, OptimalityTol = TOL^2), warmstart = true)
    # solve!(q, SCSSolver(verbose=false, suppress_warnings = true), warmstart = true)
    # solve!(q, ECOSSolver(verbose=false, abstol = 1e-6))
    solve!(q, MosekSolver(MSK_IPAR_LOG = 0, MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-3 * abs_tol))

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
    if gap < abs_tol
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
        while pernum > max_rep && lev_tol >= lev_tol_low  # println("reduce permnum!")
          lev_tol *= μ
          (c, pernum, perm) = greedy_rand(F, optsol, lev_tol)
        end
      end
      b = greedy(F, var.value[:])
      while in(b, subgradients)
        # println("b included!")
        (b, pernum1, perm) = greedy_rand(F, var.value[:], lev_tol)
      end
      var.value[:] = optsol
      while length(permset) >= pernum && lev_tol <= lev_tol_up  # println("enlarge permnum!")
        lev_tol *= λ
        (a, pernum, perm) = greedy_rand(F, optsol, lev_tol)
      end
      (a, pernum, perm) = greedy_rand(F, optsol, lev_tol)
      while in(a, subgradients)
        # println("a included!")
        if !in(perm, permset)
          push!(permset, perm)
        end
        (a, pernum) = greedy_rand(F, optsol, lev_tol)
      end
      push!(permset, perm)
    else                                           # update the upper bound
      rep = 1                                      # reset repetition count
      lev_tol = lev_tol_def                        # reset the level set tolerance
      (a, pernum, perm) = greedy_rand(F, optsol, lev_tol)
      a = greedy(F, optsol)
      permset = [sortperm(var.value[:], rev = true)]
    end

    cos_num = length(constraints)
    if cos_num > n
      prod = zeros(cos_num)
      for i = 1:cos_num
        prod[i] = dot(optsol, subgradients[i])
      end
      if cos_num == (2 * n + 3)
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

function dual_to_primal(p::Problem)
  # convert a problem in DualModel to f(x) + t
  c, A, b, dual_var_cones, var_to_ranges, vartypes, conic_constraints = conic_problem(p)
  # convert the problem into its dual form
  for i = 1:length(dual_var_cones)
    if dual_var_cones[i][1] == :Free
      dual_var_cones[i] = (:Zero, dual_var_cones[i][2])
    elseif dual_var_cones[i][1] == :Zero
      dual_var_cones[i] = (:Free, dual_var_cones[i][2])
    elseif dual_var_cones[i][1] == :ExpPrimal
      dual_var_cones[i] = (:ExpDual, dual_var_cones[i][2])
    elseif dual_var_cones[i][1] == :ExpDual
      dual_var_cones[i][1] = (:ExpPrimal, dual_var_cones[i][2])
    end
  end
  c, A, b = vec(full(b)), -A', vec(full(c))
  pos = dual_var_cones[length(dual_var_cones)].stop + 1   # the position of the new variable
  # update the paramters
  push!(c, 1)
  A = [A sparse(zeros(size(A)[1], 1))]
  push!(dual_var_cones, (:Free, pos:pos))             # the last variable is unconstrained
  dual_constr_cones = fill((:Zero, 1:size(A, 1)),1)
  return c, A, b, dual_constr_cones, dual_var_cones
end
