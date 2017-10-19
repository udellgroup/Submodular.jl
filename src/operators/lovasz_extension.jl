#############################################################################
# lovasz_extension.jl
# Handles the Lovasz extensions of set functions.
#############################################################################

export lovasz, LovaszExtAtom
export sign, monotonicity, curvature, evaluate
export ConvexLovasz

type LovaszExtAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr, Variable}
  size::Tuple{Int, Int}
  variable::Variable
  func::CombiFunc

  function LovaszExtAtom(F::CombiFunc, x::Variable)
    if length(F.setvariables) == 1
      if evaluate(F, [])[1] != 0
        error("A combinatorial function should be 0 at the empty set to derive its Lovasz extension.")
      else
        if x.size[1] == F.setvariables[1].cardinality
          children = (F, x)
          return new(:lovasz, hash(children), children, (1, 1), x, F)
        else
          error("The size of the continuous variable should be the same as the baseset of the combinatorial variable of the combinatorial function.")
        end
      end
    else
      error("A combinatorial function should be sigle-variant to obtain its Lovasz extension.")
    end
  end
end

lovasz(F::CombiFunc, x::Variable) = LovaszExtAtom(F, x)

function sign(x::LovaszExtAtom)
  return sign(x.children[1])
end

function monotonicity(x::LovaszExtAtom)
  return monotonicity(x.children[1])
end

function curvature(x::LovaszExtAtom)
  modd = modularity(x.children[1])
  if modd == SubModularity()
    return ConvexVexity()
  elseif modd == SuperModularity()
    return ConcaveVexity()
  elseif modd == Modularity()
    return AffineVexity()
  elseif modd == ConstModularity()
    return ConstVexity()
  end
end

vexity(x::LovaszExtAtom) = curvature(x)

function evaluate(F::LovaszExtAtom)
  n = F.children[2].size[1]
  y = F.children[2].value[:, 1]
  i = sortperm(y, rev = true)
  V = sort(F.children[1].setvariables[1].baseset)
  storage = zeros(n + 1)
  x = zeros(n)
  for ii = 2:n+1
    storage[ii] = evaluate(F.children[1], V[i[1:ii-1]])[1]
    x[i[ii - 1]] = storage[ii] - storage[ii - 1]
  end
  # Compatibility with Convex.jl
  val = zeros(1, 1)
  val[1] = dot(x, y)
  return val
end

function get_cv(x::LovaszExtAtom)
  variable = []
  return push!(variable, x.variable)
end

function get_sv(x::LovaszExtAtom)
  return []
end

function get_v(x::LovaszExtAtom)
  variable = []
  return push!(variable, x.variable)
end

# Convex function + Lovasz extension
type ConvexLovasz <: SCOPEModel
  convex_part::Problem
  lovasz::LovaszExtAtom
end
