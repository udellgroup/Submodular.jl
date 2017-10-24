#############################################################################
# types.jl
# Defines SubmodFunc and Sets.
# SubmodFunc subtypes of Convex.AbstractExpr, and is subtyped by atoms which
# are combinatorial functions.
# Each type which subtypes SubmodFunc must have:
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
export SubmodFunc, CombiSet, AllCombiSet, ContiSet, SCOPEModel

# Type of values
const Val = Union{Number, AbstractArray}
const ValOrNothing = Union{Val, Void}

### Combinatorial functions
abstract type SubmodFunc <: AbstractExpr end

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
    if isa(w, AbstractArray{Int}) == false
      error("The elements assigned must be integers.")
    elseif maximum(w) > maximum(var[1].baseset) || minimum(w) < minimum(var[1].baseset)
      error("The elements assigned exceeds the cardinality of the base set of the variable.")
    else
      var[1].elements = w
    end
  end
  evaluate(f)
end

function evaluate(f::AbstractExpr, w::Val...)
  var = get_v(f)
  for i = 1:length(w)
    if isa(var[1], Variable)
      var[1].value = w
    else
      if isa(w, AbstractArray{Int}) == false
        error("The elements assigned must be integers.")
      elseif maximum(w) > maximum(var[1].baseset) || minimum(w) < minimum(var[1].baseset)
        error("The elements assigned exceeds the cardinality of the base set of the variable.")
      else
        var[1].elements = w
      end
    end
  end
  evaluate(f)
end
