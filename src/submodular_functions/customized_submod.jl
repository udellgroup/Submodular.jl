#############################################################################
# customized_submod.jl
# allow users to customize submodular functions
# no check is done for the submodularity or cardinality of the function
#############################################################################

export submod
export sign, monotonicity, modularity, evaluate

type SubmodAtom <: SubmodFunc
  head::Symbol
  id_hash::UInt64
  children::Tuple{Function, SetVariable}
  size::Tuple{Int, Int}
  func::Function
  setvariables::Array{CombiSet}

  function SubmodAtom(F::Function, S::SetVariable)
    warn("Please make sure the function is submodular, the function returns 0 at the empty set, and the cardinality of the base set is equal to that of the set variable.")
    children = (F, S)
    setvariables = get_sv(S)
    return new(:custoized_submod, hash(children), children, (1, 1), F, setvariables)
  end
end

submod(F::Function, S::SetVariable) = SubmodAtom(F, S)

function sign(F::SubmodAtom)
  return NoSign()
end

function monotonicity(F::SubmodAtom)
  return NoMonotinicity()
end

function modularity(F::SubmodAtom)
  return SubModularity()
end

function evaluate(F::SubmodAtom)
  return F.func(F.children[2].elements)
end
