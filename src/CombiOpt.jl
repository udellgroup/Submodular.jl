module CombiOpt

# data structures
include("types.jl")

# combinatorial sets
include("sets/perm.jl")
include("sets/poly.jl")

# operators
include("operators/lovasz_extension.jl")

# algorithms
# include("alg/frank_wolfe.jl")

# rounding
# include("round/deterministic.jl")

end # module
