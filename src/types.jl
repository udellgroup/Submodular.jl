import Convex: Constraint

### Combinatorial sets

abstract type CombiSet end

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
