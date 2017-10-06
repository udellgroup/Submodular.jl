module CombiOpt

import MathProgBase: ConicModel, loadproblem!, optimize!, numvar

import Convex: AbstractExpr, Variable, Constraint, Constant, Solution, Problem, convert
import Convex: Vexity, ConstVexity, AffineVexity, ConvexVexity, ConcaveVexity, NotDcp
import Convex: Sign, ComplexSign, Monotonicity, Nonincreasing, Nondecreasing, NoMonotonicity
import Convex: curvature, evaluate, monotonicity, sign, vexity
import Convex: fix!, free!
import Convex: min, add_constraints!
import Convex: Sign, Positive, Negative, NoSign, ComplexSign
import Convex: +, -, *, .*, /, ./, abs, AbsAtom
import Convex: solve!
import Convex: conic_problem

export AbstractExpr, Variable, Solution, convert
export curvature, evaluate, monotonicity, sign, vexity
export +, -, *, .*, /, ./, abs
export add_constraints!
export solve!

# data structures
include("types.jl")
include("modularity.jl")
include("problems.jl")
include("variables.jl")
include("solutions.jl")

# combinatorial sets
include("combinatorial_sets/combinatorial_sets.jl")
include("combinatorial_sets/intersect.jl")
include("combinatorial_sets/setdiff.jl")
include("combinatorial_sets/union.jl")

# continuous sets
include("continuous_sets/set_constraints.jl")
include("continuous_sets/perm.jl")
include("continuous_sets/poly.jl")

# combinatorial functions
include("combinatorial_functions/card.jl")
include("combinatorial_functions/modular.jl")

# models
include("models.jl")

# operators
include("operators/affine_projection.jl")
include("operators/lovasz_extension_abs.jl")
include("operators/lovasz_extension.jl")
include("operators/greedy.jl")
include("operators/solve_dual.jl")

# algorithms
include("algos/mininum_norm_point.jl")

# prox
include("prox/prox_convex.jl")
include("prox/prox_poly.jl")

# solvers
include("solvers/chambolle_pock.jl")
include("solvers/cutting_plane.jl")
include("solvers/lp_over_assocpoly.jl")
include("solvers/proximal_level_bundle.jl")

# utilities
include("utilities/show.jl")
include("utilities/compatibilities.jl")


# import Graphs: GenericGraph, simple_graph, add_edge!, connected_components
#
# export simple_graph, add_edge!, connected_components

end # module
