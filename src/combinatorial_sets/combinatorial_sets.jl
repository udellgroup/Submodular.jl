#############################################################################
# contiuous_sets.jl
# Handles generic combinatorial sets and their operations.
#############################################################################

export GenCombiSet, get_elements, sign, evaluate

# Instance of generic Combinatorial Sets
abstract type GenCombiSet <: CombiSet end

function get_elements(x::CombiSet)
  if x.elements == nothing
  error("This set variable has not been assigned elements yet; cannot operate on unassigned set variables.")
  else
    return x.elements
  end
end

function sign(x::CombiSet)
  return x.sign
end


function evaluate(x::GenCombiSet)
  return x.value == nothing ? error("Value of the set variable is yet to be calculated") : x.value
end
