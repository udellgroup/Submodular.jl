#############################################################################
# modular.jl
# Handles modular functions of a given vector x.
#############################################################################

export modular
export sign, monotonicity, modularity, evaluate

type ModularAtom <: CombiFunc   # TODO: extend to sets other than plain set varables
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr, CombiSet}
  size::Tuple{Int, Int}
  param::AbstractExpr
  setvariables::Array{CombiSet}

  function ModularAtom(x::AbstractExpr, S::CombiSet)
    if x.size != (S.cardinality, 1)
      error("Cannot define a modular function when the vector's size is different from the (card(baseset), 1).")
    else
      children = (x, S)
      setvariables = get_sv(S)
      return new(:modular, hash(children), children, (1, 1), x, setvariables)
    end
  end
end

modular(x::Variable, S::CombiSet) = ModularAtom(x, S)

modular(x::Val, S::CombiSet) = ModularAtom(Constant(x), S)

function sign(x::ModularAtom)
  return sign(x.children[1]) * sign(x.children[2])
end

function monotonicity(x::ModularAtom)
  return (sign(x.children[2]) * Nondecreasing(), sign(x.children[1]) * Nondecreasing())
end

function modularity(x::ModularAtom)
  return Modularity()
end

function evaluate(x::ModularAtom) # TODO only allow scalar input
  summ = 0.0
  param = evaluate(x.param)
  elements = []
  elements = vcat(elements, get_elements(x.children[2]))
  elements = Set(elements)
  for i in elements
    summ += param[i]
  end
  # Compatibility with Convex.jl
  val = zeros(1, 1)
  val[1] = summ
  return val
end
