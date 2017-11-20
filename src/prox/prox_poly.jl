#############################################################################
# prox_poly.jl
# compute the prox of indicator functions of associated polyhedra
#############################################################################

export prox

function prox(p::AssocPoly, w::AbstractArray, Tol::Float64 = 1e-3; solver::String = "fujishige_wolfe")
  if solver == "card_inc_fix"
    return  card_inc_fix(p.F, w)
  elseif solver == "fujishige_wolfe"
    return fujishige_wolfe(p, w, Tol)
  elseif solver == "frank_wolfe_away"
    x = Variable(p.F.setvariables[1].cardinality)
    g = norm(x)
    prob_perm = SCOPEminimize(g, x in p)
    return frank_wolfe_away(prob_perm, verbose = false)
  elseif solver == "cutting_plane"
    x = Variable(p.F.setvariables[1].cardinality)
    q = Problem(:minimize, 0.5*norm(x)^2)
    if isa(p, BasePoly)
      f = lovasz(p.F, x)
    elseif isa(p, SymPoly)
      f = lovasz(p.F, abs(x))
    else
      error("The cutting plane method is not applicable to $(typeof(p)).")
    end
    primal_solution, primal_optval = cutting_plane(q, f)
    return -primal_solution                           # the solution we get here is the negative
  else
    error("This method is not supported.")
  end
end

prox(C::SetConstraint{AssocPoly}, w::AbstractArray, Tol::Float64 = 1e-3) = prox(p, w, Tol)
