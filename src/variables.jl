#############################################################################
# set_variables.jl
# Defines SetVariable, which is a subtype of AbstractExpr
#############################################################################

export SetVariable
export sign, evaluate, baseset, fix!, free!, elements
export get_cv, get_sv, get_v

type SetVariable <: CombiSet
  head::Symbol
  id_hash::UInt64
  elements::ValOrNothing
  baseset::AbstractArray
  cardinality::Int
  value::ValOrNothing
  sign::Sign
  sets::Array{Symbol,1}

  function SetVariable(baseset::AbstractArray, sign::Sign=NoSign(), sets::Symbol...)
    this = new(:setvariable, 0, nothing, baseset, length(baseset), nothing, sign, Symbol[sets...])
    this.id_hash = object_id(this)
    return this
  end
  SetVariable(m::Int, sign::Sign=NoSign(), sets::Symbol...) = SetVariable(collect(1:m), sign, sets...)
  SetVariable(sign::Sign, sets::Symbol...) = SetVariable([], sign, sets...)
  SetVariable(sets::Symbol...) = SetVariable([], NoSign(), sets...)
  SetVariable(baseset::AbstractArray, sets::Symbol...) = SetVariable(baseset, NoSign(), sets...)
  SetVariable(m::Int, sets::Symbol...) = SetVariable(collect(1:m), sets...)
end

id_to_variables = Dict{UInt64, SetVariable}()

function evaluate(x::SetVariable)
  return x.value == nothing ? error("Value of the set variable is yet to be calculated") : x.value
end

function evaluate(x::SetVariable, w::Val)
  x.elements = w
end

function sign(x::SetVariable)
  return x.sign
end

# fix set variables to hold them at their current value, and free them afterwards
function fix!(x::SetVariable)
  x.value == nothing && error("This set variable has not been assigned elements yet; cannot fix value to nothing!")
  push!(x.sets, :fixed)
  x
end

function fix!(x::SetVariable, v)
  # TODO: check value inclusion
  x.value = v
  fix!(x)
end

function free!(x::SetVariable)
  # TODO this won't work if :fixed appears other than at the end of x.sets
  x.sets[end] == :fixed && pop!(x.sets)
  x
end

### Output the variables of an expression
function get_cv(x::AbstractExpr)
  if isa(x, Variable) || (typeof(x) == Constant) || (typeof(x) == SetVariable)
    if isa(x, Variable)
      cv = []
      push!(cv, x)
      return cv
    else
      return []
    end
  else
    cv = []
    for i = 1:length(x.children)
      cvi = get_cv(x.children[i])
      cv = vcat(cv, cvi)
    end
    cvSet = Set(cv)
    cv = collect(cvSet)
    return cv
  end
end

function get_cv(x::Val)
  return []
end

### Output the set variables of an expression
function get_sv(x::AbstractExpr)
  if isa(x) == Variable || (typeof(x) == Constant) || (typeof(x) == SetVariable)
    if isa(x) == SetVariable
      sv = []
      push!(sv, x)
      return sv
    else
      return []
    end
  else
    sv = []
    for i = 1:length(x.children)
      svi = get_sv(x.children[i])
      sv = vcat(sv, svi)
    end
    svSet = Set(sv)
    sv = collect(svSet)
    return sv
  end
end

function get_sv(f::CombiFunc)
  return f.setvariables
end

function get_sv(x::Val)
  return []
end

### Output all the variables of an expression
function get_v(x::AbstractExpr)
  if isa(x, Variable) || (typeof(x) == Constant) || (typeof(x) == SetVariable)
    if isa(x, Variable) || typeof(x) == SetVariable
      v = []
      push!(v, x)
      return v
    else
      return []
    end
  else
    v = []
    for i = 1:length(x.children)
      vi = get_v(x.children[i])
      v = vcat(v, vi)
    end
    vSet = Set(v)
    v = collect(vSet)
    return v
  end
end

function get_v(f::CombiFunc)
  return f.setvariables
end

function get_v(x::Val)
  return []
end
