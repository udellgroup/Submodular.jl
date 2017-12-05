module Submodular

import Convex: AbstractExpr, Variable, Constraint, Constant, Solution, Problem, convert
import Convex: Vexity, ConstVexity, AffineVexity, ConvexVexity, ConcaveVexity, NotDcp
import Convex: Sign, ComplexSign, Monotonicity, Nonincreasing, Nondecreasing, ConstMonotonicity, NoMonotonicity
import Convex: curvature, evaluate, monotonicity, sign, vexity
import Convex: fix!, free!
import Convex: min, add_constraints!
import Convex: Sign, Positive, Negative, NoSign, ComplexSign
import Convex: +, -, *, .*, /, ./, abs, AbsAtom, vecdot
import Convex: minimize, maximize, solve!
import Convex: conic_problem

import ForwardDiff: gradient

import LightGraphs: AbstractGraph, AbstractEdge, Graph, nv, ne, edges, Edge, add_edge!, weights, src, dst, induced_subgraph, connected_components, modularity

import MathProgBase: ConicModel, loadproblem!, optimize!, numvar, AbstractMathProgSolver

export AbstractExpr, Variable, Solution, convert
export curvature, evaluate, monotonicity, sign, vexity
export +, -, *, .*, /, ./, abs
export add_constraints!
export solve!

export conic_problem

# data structures
include("types.jl")
include("modularity.jl")
include("graphs.jl")
include("problems.jl")
include("variables.jl")


# combinatorial sets
include("combinatorial_sets/combinatorial_sets.jl")
include("combinatorial_sets/intersect.jl")
include("combinatorial_sets/setdiff.jl")
include("combinatorial_sets/union.jl")

# continuous sets
include("continuous_sets/set_constraints.jl")
include("continuous_sets/perm.jl")
include("continuous_sets/poly.jl")

# submodular functions
include("submodular_functions/card_based.jl")
include("submodular_functions/cut.jl")
include("submodular_functions/log_determinant.jl")
include("submodular_functions/rank_of_graph_matroid.jl")
include("submodular_functions/modular.jl")
include("submodular_functions/customized_submod.jl")

# operators
include("operators/affine_projection.jl")
include("operators/dual_to_primal.jl")
include("operators/extreme_point.jl")
include("operators/lovasz_extension_abs.jl")
include("operators/lovasz_extension.jl")
include("operators/gradient.jl")
include("operators/greedy.jl")

# prox
include("prox/prox_convex.jl")
include("prox/prox_poly.jl")

# models
include("models.jl")

# solutions
include("solutions.jl")

# algorithms
include("algorithms/card_inc_fix.jl")
include("algorithms/chambolle_pock.jl")
include("algorithms/cutting_plane.jl")
include("algorithms/frank_wolfe_away.jl")
include("algorithms/fujishige_wolfe.jl")
include("algorithms/lp_over_assocpoly.jl")
include("algorithms/L_BFGS.jl")
# include("algorithms/proximal_level_bundle.jl")

# utilities
include("utilities/show.jl")

end # module
