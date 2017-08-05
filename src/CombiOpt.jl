module CombiOpt

# data structures
include("types.jl")

# combinatorial sets
include("sets/perm.jl")
include("sets/poly.jl")

# operators
include("operators/aff_projection.jl")
include("operators/lovasz_extension.jl")
include("operators/greedy.jl")

# algorithms
include("algos/minimum_norm.jl")

# prox
include("prox/prox_poly.jl")

# rounding
# include("round/deterministic.jl")

end # module
