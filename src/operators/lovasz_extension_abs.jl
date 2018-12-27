#############################################################################
# lovasz_extension.jl
# Handles the Lovasz extensions of set functions.
#############################################################################

export lovasz, LovaszExtAbsAtom
export sign, monotonicity, curvature, evaluate
export ConvexLovaszAbs

mutable struct LovaszExtAbsAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr, Variable}
  size::Tuple{Int, Int}
  variable::Variable
  func::SubmodFunc

  function LovaszExtAbsAtom(F::SubmodFunc, x::AbsAtom)
    if length(F.setvariables) == 1
      if length(x.children) == 1 && isa(x.children[1], Variable)
        if evaluate(F, [])[1] != 0
          error("A combinatorial function should be 0 at the empty set to derive its Lovasz extension.")
        else
          var = get_v(x)
          if length(var) != 1
            error("Functions with other than one variable are not supported.")
          else
            if var[1].size[1] == F.setvariables[1].cardinality
              children = (F, var[1])
              return new(:lovaszabs, hash(children), children, (1, 1), var[1], F)
            else
              error("The size of the continuous variable should be the same as the baseset of the combinatorial variable of the combinatorial function.")
            end
          end
        end
      else
        error("Only able to define Lovasz extesions on absolute values of a variable.")
      end
    else
      error("A combinatorial function should be sigle-variant to obtain its Lovasz extension.")
    end
  end
end

lovasz(F::SubmodFunc, x::AbsAtom) = LovaszExtAbsAtom(F, x)

function sign(F::LovaszExtAbsAtom)
  return sign(F.children[1])
end

function monotonicity(F::LovaszExtAbsAtom)
  return monotonicity(F.children[1])
end

function curvature(F::LovaszExtAbsAtom)
  modd = modularity(F.children[1])
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

vexity(F::LovaszExtAbsAtom) = curvature(F)

function evaluate(f::LovaszExtAbsAtom)
  n = f.children[2].size[1]
  y = abs.(f.children[2].value[:])
  i = sortperm(y, rev = true)
  V = sort(f.children[1].setvariables[1].baseset)
  storage = zeros(n + 1)
  x = zeros(n)
  for ii = 2:n+1
    storage[ii] = evaluate(f.children[1], V[i[1:ii-1]])[1]
    x[i[ii - 1]] = storage[ii] - storage[ii - 1]
  end
  # Compatibility with Convex.jl
  val = zeros(1, 1)
  val[1] = dot(x, y)
  return val
end

function get_cv(F::LovaszExtAbsAtom)
  variable = []
  return push!(variable, F.variable)
end

function get_sv(F::LovaszExtAbsAtom)
  return []
end

function get_v(F::LovaszExtAbsAtom)
  variable = []
  return push!(variable, F.variable)
end
