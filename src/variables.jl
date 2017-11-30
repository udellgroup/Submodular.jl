#############################################################################
# set_variables.jl
# Defines SetVariable, which is a subtype of AbstractExpr
#############################################################################

export SetVariable
export fix!, free!, elements
export get_cv, get_sv, get_v

type SetVariable <: CombiSet
  head::Symbol
  id_hash::UInt64
  elements::ValOrNothing
  baseset::AbstractArray
  cardinality::Int

  function SetVariable(baseset::AbstractArray)
    this = new(:setvariable, 0, nothing, baseset, length(baseset))
    this.id_hash = object_id(this)
    return this
  end
  SetVariable(m::Int) = SetVariable(collect(1:m))
end

id_to_variables = Dict{UInt64, SetVariable}()

# fix set variables to hold them at their current value, and free them afterwards
function fix!(x::SetVariable)
  x.elements == nothing && error("This set variable has not been assigned elements yet; cannot fix value to nothing!")
  x
end

function fix!(x::SetVariable, v)
  # TODO: check value inclusion
  x.elements = v
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
      cv = Variable[]
      push!(cv, x)
      return cv
    else
      return Variable[]
    end
  else
    cv = Variable[]
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
  return Variable[]
end

### Output the set variables of an expression
function get_sv(x::AbstractExpr)
  if isa(x, Variable) || (typeof(x) == Constant) || (typeof(x) == SetVariable)
    if isa(x, SetVariable)
      sv = SetVariable[]
      push!(sv, x)
      return sv
    else
      return SetVariable[]
    end
  else
    sv = SetVariable[]
    for i = 1:length(x.children)
      svi = get_sv(x.children[i])
      sv = vcat(sv, svi)
    end
    svSet = Set(sv)
    sv = collect(svSet)
    return sv
  end
end

function get_sv(f::SubmodFunc)
  return f.setvariables
end

function get_sv(x::Val)
  return []
end

### Output all the variables of an expression
function get_v(x::AbstractExpr)
  if isa(x, Variable) || (typeof(x) == Constant) || (typeof(x) == SetVariable)
    if isa(x, Variable) || typeof(x) == SetVariable
      v = AbstractExpr[]
      push!(v, x)
      return v
    else
      return []
    end
  else
    v = AbstractExpr[]
    for i = 1:length(x.children)
      vi = get_v(x.children[i])
      v = vcat(v, vi)
    end
    vSet = Set(v)
    v = collect(vSet)
    return v
  end
end

function get_v(f::SubmodFunc)
  return f.setvariables
end

function get_v(x::Val)
  return []
end
