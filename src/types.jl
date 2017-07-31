import Convex: Constraint
import Convex: AbstractExpr

### Combinatorial sets

abstract type CombiSet end

### Combinatorial functions
abstract type CombiFunc <: AbstractExpr end

# the type of generic combinatorial fucntions
type GenCombi <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr}
  size::Tuple{Int, Int}
end

### Combinatorial problems

# a CombiProblem represents the combinatorial problem
#   minimize objective(x)
#   st       x \in s, s in combi_constraints
#            x \in c, c in convex_constraints


# type CombiProblem{F}
#   objective::F # F could be DiffableFunction, AbstractExpr, ...
#   combi_constraints::Array{CombiSet}
#   convex_constraints::Array{Constraint}
# end
