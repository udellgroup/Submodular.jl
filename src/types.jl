import MathProgBase
import Convex: Constraint, AbstractExpr

export CombiFunc, GenCombi, CombiSet, Problem, Solution
export Float64OrNothing

### Combinatorial functions
abstract type CombiFunc <: AbstractExpr end

# the type of generic combinatorial fucntions
type GenCombi <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr}
  size::Tuple{Int, Int}
end

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

const Float64OrNothing = Union{Float64, Void}

# TODO: Cleanup
type Solution{T<:Number}
  primal::Array{T, 1}
  dual::Array{T, 1}
  status::Symbol
  optval::T
  has_dual::Bool
end

Solution{T}(x::Array{T, 1}, status::Symbol, optval::T) = Solution(x, T[], status, optval, false)
Solution{T}(x::Array{T, 1}, y::Array{T, 1}, status::Symbol, optval::T) = Solution(x, y, status, optval, true)

type Problem
  head::Symbol
  objective::AbstractExpr
  combi_constraints::Array{CombiSet}
  convex_constraints::Array{Constraint}
  status::Symbol
  optval::Float64OrNothing
  model::MathProgBase.AbstractConicModel
  solution::Solution

  function Problem(head::Symbol, objective::AbstractExpr,
                   model::MathProgBase.AbstractConicModel, constraints::Array=Constraint[])
    if sign(objective)== Convex.ComplexSign()
      error("Objective can not be a complex expression")
    else
      return new(head, objective, constraints, "not yet solved", nothing, model)
    end
  end
end
