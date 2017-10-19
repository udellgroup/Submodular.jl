import JuMP: DiffableFunction

function frank_wolfe(p::CombiProblem{DiffableFunction})
  # check only constraint is one combinatorial constraint
  @assert length(p.convex_constraints==0) && length(p.combi_constraints==1)
  combiset = p.combi_constraints[1]
  n = combiset.dim
  x = zeros(n)
  for k=1:100 # maxiters
    g = grad(p.objective, x)
    h = fenchel(combiset, g)
    alpha = 1/k # stepsize
    x = (1-alpha)*x + alpha*h
  end
  return x
end
