#############################################################################
# compatibilities.jl
# Deals with SCOPE.jl's compatibilities with Convex.jl.
#############################################################################

export zer0

# Outputs of Convex.jl are AbstractMatrix instead of AbstractArray
const zer0 = zeros(1, 1)
