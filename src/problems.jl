#############################################################################
# type SCOPEProblem{M}
#   head::Symbol                      # :minimize or :maximize
#   objective::AbstractExpr
#   status::Symbol                    # the status of the problem
#   optval::Float64OrNothing
#   model::T                          # a subtype of PrimalModel or DualModel, see models.jl
#   solution::Solution
#############################################################################

export SCOPEProblem, SCOPEminimize, SCOPEmaximize, satisfy, add_constraints!

const Float64OrNothing = Union{Float64, Nothing}

mutable struct SCOPEProblem{T <: SCOPEModel}
  head::Symbol
  objective::AbstractExpr
  constraints::AbstractArray
  status::String
  optval::Float64OrNothing
  model::T
  # solution::Solution
end

function SCOPEProblem(head::Symbol, objective::AbstractExpr, model::SCOPEModel,
                 constraints::Array=Constraint[])
  if sign(objective) == ComplexSign()
    error("Objective can not be a complex expression")
  else
    return SCOPEProblem(head, objective, constraints, "not yet solved", nothing, model)
  end
end

# constructor if model is not specified
function SCOPEProblem(head::Symbol, objective::AbstractExpr,
                 constraints::AbstractArray=Constraint[],
                 model::SCOPEModel = get_model(objective, constraints))
  SCOPEProblem(head, objective, model, constraints)
end

SCOPEProblem(head::Symbol, objective::AbstractExpr, constraints::Constraint...) =
  SCOPEProblem(head, objective, [constraints...])

# Allow users to simply type minimize
SCOPEminimize(objective::AbstractExpr, constraints::Constraint...) =
  SCOPEProblem(:minimize, objective, collect(constraints))
SCOPEminimize(objective::AbstractExpr, constraints::Array{<:Constraint}=Constraint[]) =
  SCOPEProblem(:minimize, objective, constraints)
SCOPEminimize(objective::Values, constraints::Constraint...) =
  minimize(convert(AbstractExpr, objective), collect(constraints))
SCOPEminimize(objective::AbstractExpr, constraints::Array{<:Constraint}=Constraint[]) =
  minimize(convert(AbstractExpr, objective), constraints)

# Allow users to simply type maximize
SCOPEmaximize(objective::AbstractExpr, constraints::Constraint...) =
  SCOPEProblem(:maximize, objective, collect(constraints))
SCOPEmaximize(objective::AbstractExpr, constraints::Array{<:Constraint}=Constraint[]) =
  SCOPEProblem(:maximize, objective, constraints)
SCOPEmaximize(objective::Values, constraints::Constraint...) =
  maximize(convert(AbstractExpr, objective), collect(constraints))
SCOPEmaximize(objective::AbstractExpr, constraints::Array{<:Constraint}=Constraint[]) =
  maximize(convert(AbstractExpr, objective), constraints)

# # Allow users to simply type satisfy (if there is no objective)
# satisfy(constraints::Constraint...) = SCOPEProblem(:minimize, Constant(0), [constraints...])
# satisfy{T<:Constraint}(constraints::Array{T}=Constraint[]) =
#   SCOPEProblem(:minimize, Constant(0), constraints)
# satisfy(constraint::Constraint) = satisfy([constraint])

# +(constraints, constraints) is defined in constraints.jl
add_constraints!(objective::AbstractExpr, constraints::Array{<:Constraint}=Constraint[]) = +(p.constraints, constraints)
add_constraints!(p::SCOPEProblem, constraint::Constraint) = add_constraints!(p, [constraint])
add_constraint! = add_constraints!
