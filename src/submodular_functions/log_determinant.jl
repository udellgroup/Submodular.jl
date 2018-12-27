#############################################################################
# log_determinant.jl
# given V = {1, 2, ..., n}, M ∈ ℜ^{n×n}, F(S) = logdet(M_S)
#############################################################################

import Base.logdet

export logdet
export sign, monotonicity, modularity, evaluate

mutable struct LogDeterminantAtom <: SubmodFunc
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractMatrix, CombiSet}
  size::Tuple{Int, Int}
  matrix::AbstractMatrix
  setvariables::Array{CombiSet}

  function LogDeterminantAtom(matrix::AbstractMatrix, S::CombiSet)
    if size(matrix)[1] != S.cardinality
      error("Cannot define a log determinant function when the number of volumns of the matrix is different from the size of the set variable.")
    else
      children = (matrix, S)
      setvariables = get_sv(S)
      return new(:logdet, hash(children), children, (1, 1), matrix, setvariables)
    end
  end
end

logdet(matrix::AbstractMatrix, S::CombiSet) = LogDeterminantAtom(matrix, S)

function sign(F::LogDeterminantAtom)
  return NoSign()
end

function monotonicity(F::LogDeterminantAtom)
  return (NoMonotonicity(), )
end

function modularity(F::LogDeterminantAtom)
  if isposdef(F.matrix)
    return SubModularity()
  elseif isposdef(-F.matrix)
    return SuperModularity()
  else
    return NotDetermined()
  end
end

function evaluate(F::LogDeterminantAtom)
  set = get_elements(F.children[2])
  if length(set) == 0
    return 0
  end
  matrix = F.matrix[]
  return logdet(F.matrix[set, set])
end
