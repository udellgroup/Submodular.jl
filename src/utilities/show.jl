import Base.show
export show

# A set variable, for example, Variable([1, 3]), will be displayed as:
# SetVariable of
# baseset: [1, 3]
# sign: NoSign()
# modularity: ConstModularity()
function show(io::IO, x::SetVariable)
  print(io, """SetVariable of
    baseset: $(x.baseset)""")
  if x.elements != nothing
    print(io, "\nelements: $(x.elements)")
  end
end

# A combinatorial set, for example, GenCombiSet([1, 3]), will be displayed as:
# CombiSet of
# baseset: [1, 3]
# sign: NoSign()
# modularity: ConstModularity()
function show(io::IO, x::CombiSet)
  print(io, """CombiSet of
    baseset: $(x.baseset)""")
  if x.elements != nothing
    print(io, "\nelements: $(x.elements)")
  end
end

# An combinatorial function, for example, card(x)^2, will be displayed as:
# SubmodFunc with
# head: card
# size: (1, 1)
# sign: Positive()
# modularity: SuperModularity()
function show(io::IO, f::SubmodFunc)
  print(io, """SubmodFunc with
    head: $(f.head)
    size: ($(f.size[1]), $(f.size[2]))
    sign: $(sign(f))
    modularity: $(modularity(f))
    """)
end

# An polyhedron asssociated with f, for example, subpoly(f), will be dispayed as:
# AssocPoly with
# head: card
# baseset: [1, 2]
# function: f
function show(io::IO, p::AssocPoly)
  print(io, """AssocPoly with
    head: $(p.head)
    baseset: $(p.V)
    associated to $(p.F)""")
end

# A set constraint, for example, in(x, p), will be displayed as:
# Constraint:
# set constraint
# lhs: ...
# rhs: ...
function show(io::IO, c::SetConstraint)
  print(io, """Constraint:
    $(c.head) constraint
    set: $(c.rhs)""")
end
