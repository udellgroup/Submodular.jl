in general the code structure is beautiful: very well organized and extremely clear

in frank wolfe solver:

* implement linesearch
* implement stopping condition
* allow users to set parameters (verbosity, maxiters, abstol, reltol)

in cutting plane solver:

* parameters should be keyword arguments to function, not global (or local) variables
* you'll want to go through when you're done debugging and make sure the names of variables make sense and are all well commented. eg "haha" is not a great variable name :)
* check your algorithm for places where you might be using more memory than you need. for example, can you reuse the same arrays at each iteration for your working variables, rather than creating new arrays each time?

general code structure:

* why call values "Vall" and not "Val" or "Value"?
* how to you handle minimization vs maximization? i don't see any tooling for this in, eg, the frank wolfe solver

julia tricks:

* to check types, use the function "isa": for example,
    `isa(var[1], Variable)` instead of `typeof(var[1]) == Variable`
