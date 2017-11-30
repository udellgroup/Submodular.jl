draft of paper is [on sharelatex](https://www.sharelatex.com/project/59b977871683de3739c3403a)

in frank wolfe solver:

* implement linesearch
* implement stopping condition
* allow users to set parameters (verbosity, maxiters, abstol, reltol)

in cutting plane solver:

* ~~parameters should be keyword arguments to function, not global (or local) variables~~
* ~~you'll want to go through when you're done debugging and make sure the names of variables make sense and are all well commented. eg "haha" is not a great variable name :)~~
* check your algorithm for places where you might be using more memory than you need. for example, can you reuse the same arrays at each iteration for your working variables, rather than creating new arrays each time?

general code structure:

* make documentation
* conversion of the dual form to the primal for the cutting plane (bundle) method
* add the BFGS method
* change name of the package
* how to you handle minimization vs maximization? i don't see any tooling for this in, eg, the frank wolfe solver
* ~~i cannot figure out what the right syntax is for forming a combinatorial problem from the code! how do i form any associated polyhedron for a submodular function? perhaps you forgot to check in that code? or perhaps i just need to see more examples...!~~
* ~~why call values "Vall" and not "Val" or "Value"?~~
* ~~fix fw with away~~
* ~~add a blackbox atom for submodular functions~~

julia tricks:

* make sure to add all packages used by this package to the `REQUIRE`
* ~~to check types, use the function "isa": for example,
    `isa(var[1], Variable)` instead of `typeof(var[1]) == Variable`~~

syntax:

* change notations, F() for submodular functions and f() for lovasz extensions.
* ~~use `minimize(objective)` rather than `min(objective)` to form a problem (mimics Convex.jl syntax, and follows the math convention that distinguishes between minimize (which forms a problem), and min (which returns a value))~~
