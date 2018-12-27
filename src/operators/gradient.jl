#############################################################################
# gradient.jl
# calculates the gradients of functions
#############################################################################

export grad
export sign, curvature, monotonicity, evaluate

mutable struct GradientAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr}
  size::Tuple{Int, Int}
  variable::Variable
  func::Function

  function GradientAtom(x::AbstractExpr)
    var = get_cv(x)
    if x.size != (1, 1) || sign(x) == ComplexSign()
      error("A function has to be real-valued to calculate its gradient.")
    elseif length(var) > 1
      error("Only able to calculate gradients of functions with single variable.")
    elseif var[1].size[2] != 1
      error("Only able to calculate gradients of functions with variable's size being (n, 1).")
    end
    children = (x,)
    f(z) = evaluate(x, z)[1]
    return new(:grad, hash(children), children, var[1].size, var[1], f)
  end
end

grad(x::AbstractExpr) = GradientAtom(x)

function sign(x::GradientAtom)
  if monotonicity(x.children[1]) == (Nonincreasing(), )
    return Positive()
  elseif monotonicity(x.children[1]) == (Nondecreasing(), )
    return Negative()
  elseif monotonicity(x.children[1]) == (ConstMonotonicity(), )
    return Positive()
  else
    return NoSign()
  end
end

function curvature(x::GradientAtom)
  return NotDcp()
end

function monotonicity(x::GradientAtom)
  if curvature(x.children[1]) == ConvexVexity()
    return (Nondecreasing(), )
  elseif curvature(x.children[1]) == ConcaveVexity()
    return (Nonincreasing(), )
  elseif curvature(x.children[1]) == AffineVexity
    return (ConstMonotonicity(), )
  else
    return (NoMonotonicity(), )
  end
end

function evaluate(x::GradientAtom)
  var = copy(x.variable.value)
  val = gradient(x.func, var)
  x.variable.value = var
  return val
end
