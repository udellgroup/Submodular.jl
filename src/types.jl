#############################################################################
# types.jl
# Defines CombiFunc and Sets.
# CombiFunc subtypes of Convex.AbstractExpr, and is subtyped by atoms which
# are combinatorial functions.
# Each type which subtypes CombiFunc must have:
#
## head::Symbol                  -- a symbol such as :vecnorm, :+ etc
## children::(AbstractExpr,)     -- the expressions on which the current expression
##                               -- is operated
## id_hash::UInt64               -- identifier hash, can be a hash of children
##                                  or a unique identifier of the object
## size::(Int, Int)              -- size of the resulting expression.
## setvariables::SetVariable      -- the set variables of this combinatorial function.
#
#
# In addition, each atom must implement the following functions:
## sign: Returns the sign of the result of the expression;
## monotonicity: The monotonicity of the arguments with respect to the function,
##      i.e if the argument is nondecreasing, will the function be nonincreasing
##      or nondecreasing? eg. negate(x) will have Nonincreasing monotonicity;
## exprmodularity: Returns the modularity of the result of the expression;
## evaluate: Evaluates the value of the expression, assuming the problem has been
##           solved.
#
# Sets is subtyped by combinatorial sets.
#############################################################################

export Val, ValOrNothing
export CombiFunc, CombiSet, AllCombiSet, ContiSet, SCOPEModel

# Type of values
const Val = Union{Number, AbstractArray}
const ValOrNothing = Union{Val, Void}

### Combinatorial functions
abstract type CombiFunc <: AbstractExpr end

### Combinatorial sets
abstract type CombiSet <: AbstractExpr end
### All combinatorial Sets
const AllCombiSet = Union{AbstractArray, CombiSet}

### Sets on continuous variables
abstract type ContiSet <: AbstractExpr end

### Problem models
abstract type SCOPEModel end

# only works for expressions with one variable
function evaluate(f::AbstractExpr, w::AbstractArray)
  var = get_v(f)
  if isa(var[1], Variable)
    var[1].value = w
  else
    var[1].elements = w
  end
  evaluate(f)
end

function evaluate(f::AbstractExpr, w::Val...)
  var = get_v(f)
  for i = 1:length(w)
    if isa(var[i], Variable)
      var[i].value = w[i]
    else
      var[i].elements = w[i]
    end
  end
  evaluate(f)
end
