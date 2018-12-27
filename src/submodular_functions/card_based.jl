#############################################################################
# card.jl
# Handles functions of the form f(A) = g(length(A)), where g is a concave function
#############################################################################

export card, compose, CardBasedAtom
export in
export sign, monotonicity, modularity, evaluate

mutable struct CardBasedAtom{T} <: SubmodFunc
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr, SetVariable}
  size::Tuple{Int, Int}
  setvariables::Array{SetVariable}
  func::T
end

function CardBasedAtom(F::Function, S::CombiSet)
  z = Variable(1)
  children = (F(z), S)                           # create an expression
  if F(0) != 0
    error("A function has to return 0 at 0 in order to derive a cardinality-based function.")
  end
  setvariables = get_sv(S)
  return(CardBasedAtom(:card, hash(children), children, (1, 1), setvariables, F))
end

function CardBasedAtom(w::Array{Float64}, S::CombiSet)
  children = (Constant(w), S)                           # create an expression
  setvariables = get_sv(S)
  if length(setvariables) != 1
    error("The expression should have exactly one set variable.")
  elseif size(w) != (setvariables[1].cardinality, )
    error("The sizes of the set variable and the array mismatch.")
  end
  w₁ = copy(w)
  for i = 2:length(w)
    w₁[i] = w[i] - w[i-1]
  end
  return(CardBasedAtom(:card, hash(children), children, (1, 1), setvariables, w₁))
end

card(f::Function, S::CombiSet) = CardBasedAtom(f, S)

card(w::Array{Float64}, S::CombiSet) = CardBasedAtom(w, S)

card(w::Array{Int}, S::CombiSet) = CardBasedAtom(Float64.(w), S)

card(S::CombiSet) = card(x -> x, S)

function sign(F::CardBasedAtom)
  return sign(F.children[1])
end

function sign(F::CardBasedAtom{Array{Float64, 1}})
  if all(x -> x>=0, F.func)
    return Positive()
  elseif all(x -> x<=0, F.func)
    return Negative()
  else
    return NoSign()
  end
end

function monotonicity(F::CardBasedAtom)
  return (monotonicity(F.children[1]), )
end

function monotonicity(F::CardBasedAtom{Array{Float64, 1}})
  if all(x -> x>=0, F.func)
    return Nondecreasing()
  elseif all(x -> x<=0, F.func)
    return Nonincreasing()
  elseif all(x -> x==0, F.func)
    return ConstMonotonicity()
  else
    return NoMonotonicity()
  end
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

function modularity(F::CardBasedAtom{Array{Float64, 1}})
  w₁ = sort(F.func, rev=true)
  if w₁ == F.func
    return SubModularity()
  else
    w₂ = sort(F.func)
    if w₂ == F.func
      return SuperModularity()
    elseif all(x -> x==0, F.func)
      return ConstModularity()
    else
      return NotDetermined()
    end
  end
end

function evaluate(F::CardBasedAtom)
  var = F.setvariables[1]
  var1 = get_v(F.children[1])
  var1[1].value = length(Set(get_elements(var)))
  evaluate(F.children[1])
end

function evaluate(F::CardBasedAtom{Array{Float64, 1}})
  return sum(F.func[1:length(get_elements(F.children[2]))])
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
    y[i + 1] = evaluate(p.F, S[push!(ind, ordering[i])])[1]
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
