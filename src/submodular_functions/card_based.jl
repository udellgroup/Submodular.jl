#############################################################################
# card.jl
# Handles functions of the form f(A) = g(length(A)), where g is a concave function
#############################################################################

export card, compose
export in
export sign, monotonicity, modularity, evaluate

type CardBasedAtom <: SubmodFunc
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr}
  size::Tuple{Int, Int}
  setvariables::Array{SetVariable}

  function CardBasedAtom(F::Function, S::CombiSet)
    z = Variable(1)
    children = (F(z),)                           # create an expression
    if F(0) != 0
      error("A function has to return 0 at 0 in order to derive a cardinality-based function.")
    end
    setvariables = get_sv(S)
    return new(:card, hash(children), children, (1, 1), setvariables)
  end
end

card(f::Function, S::CombiSet) = CardBasedAtom(f, S)

card(S::CombiSet) = card(x -> x, S)

function sign(F::CardBasedAtom)
  return sign(F.children[1])
end

function monotonicity(F::CardBasedAtom)
  return (monotonicity(F.children[1]), )
end

function modularity(F::CardBasedAtom)
  if vexity(F.children[1]) == ConcaveVexity()
    return SubModularity()
  elseif vexity(F.children[1]) == ConvexVexity()
    return SuperModularity()
  elseif vexity(F.children[1]) == AffineVexity()
    return Modularity()
  elseif vexity(F.children[1]) == ConstVexity()
    return ConstModularity()
  else
    return NotDetermined()
  end
end

function evaluate(F::CardBasedAtom)
  var = F.setvariables[1]
  var1 = get_v(F.children[1])
  var1[1].value = length(Set(get_elements(var)))
  evaluate(F.children[1])
end

function in(w::AbstractArray, p::SubmodPoly{CardBasedAtom})
  n = length(w)
  @assert length(p.V) == n
  checker = true
  S = sort(p.V)
  ordering = sortperm(w, rev = true)
  h = 0
  ind = []
  y = zeros(n + 1)
  for i = 1:length(w)
    y[i + 1] = evaluate(p.f, S[push!(ind, ordering[i])])[1]
    h += - y[i] - w[i] + y[i + 1]
    if h < 0
      checker = false
      break
    end
  end
  return checker
end

function compose(f::Function, g::CardBasedAtom)
  return card(f, g.setvariables[1])
end
