module CombiOpt

# data structures
include("types.jl")

# combinatorial sets
include(set/perm.jl)

# algorithms
include("alg/frank_wolfe.jl")

# rounding
# include("round/deterministic.jl")

end # module
