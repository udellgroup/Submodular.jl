#############################################################################
# a CombiProblem represents the combinatorial problem
#   minimize objective(x)
#   st       x \in s, s in combi_constraints
#            x \in c, c in convex_constraints

# type CombiProblem{F}
#   objective::F # F could be DiffableFunction, AbstractExpr, ...
#   combi_constraints::Array{Sets}
#   convex_constraints::Array{Constraint}
# end
#############################################################################

export SCOPEPrimal, SCOPEminimize, SCOPEmaximize, satisfy, add_constraints!

const Float64OrNothing = Union{Float64,Void}

type SCOPEPrimal
  head::Symbol
  objective::AbstractExpr
  constraints::Array{Constraint}
  status::Symbol
  optval::Float64OrNothing
  model::SCOPEModel
  solution::Solution

  function SCOPEPrimal(head::Symbol, objective::AbstractExpr, model::SCOPEModel,
                   constraints::Array=Constraint[])
    if sign(objective) == ComplexSign()
      error("Objective can not be a complex expression")
    else
      return new(head, objective, constraints, "not yet solved", nothing, model)
    end
  end
end

# constructor if model is not specified
function SCOPEPrimal(head::Symbol, objective::AbstractExpr,
                 constraints::AbstractArray=Constraint[],
                 model::SCOPEModel = get_model(objective, constraints))
  SCOPEPrimal(head, objective, model, constraints)
end

SCOPEPrimal(head::Symbol, objective::AbstractExpr, constraints::Constraint...) =
  SCOPEPrimal(head, objective, [constraints...])

# Allow users to simply type minimize
SCOPEminimize(objective::AbstractExpr, constraints::Constraint...) =
  SCOPEPrimal(:minimize, objective, collect(constraints))
SCOPEminimize{T<:Constraint}(objective::AbstractExpr, constraints::Array{T}=Constraint[]) =
  SCOPEPrimal(:minimize, objective, constraints)
SCOPEminimize(objective::Val, constraints::Constraint...) =
  minimize(convert(AbstractExpr, objective), collect(constraints))
SCOPEminimize{T<:Constraint}(objective::Val, constraints::Array{T}=Constraint[]) =
  minimize(convert(AbstractExpr, objective), constraints)

# Allow users to simply type maximize
SCOPEmaximize(objective::AbstractExpr, constraints::Constraint...) =
  SCOPEPrimal(:maximize, objective, collect(constraints))
SCOPEmaximize{T<:Constraint}(objective::AbstractExpr, constraints::Array{T}=Constraint[]) =
  SCOPEPrimal(:maximize, objective, constraints)
SCOPEmaximize(objective::Val, constraints::Constraint...) =
  maximize(convert(AbstractExpr, objective), collect(constraints))
SCOPEmaximize{T<:Constraint}(objective::Val, constraints::Array{T}=Constraint[]) =
  maximize(convert(AbstractExpr, objective), constraints)

# # Allow users to simply type satisfy (if there is no objective)
# satisfy(constraints::Constraint...) = SCOPEPrimal(:minimize, Constant(0), [constraints...])
# satisfy{T<:Constraint}(constraints::Array{T}=Constraint[]) =
#   SCOPEPrimal(:minimize, Constant(0), constraints)
# satisfy(constraint::Constraint) = satisfy([constraint])

# +(constraints, constraints) is defined in constraints.jl
add_constraints!{T<:Constraint}(p::SCOPEPrimal, constraints::Array{T}) = +(p.constraints, constraints)
add_constraints!(p::SCOPEPrimal, constraint::Constraint) = add_constraints!(p, [constraint])
add_constraint! = add_constraints!
