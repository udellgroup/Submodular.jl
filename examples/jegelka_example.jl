using HDF5, JLD
using Distributions

#constructing a random bipartite graph (U = N, V= N) with degree varying from
#~0.15 N to ~0.65 N (effectively we disregard overlapping edges). Higher degree will
#result is a greater overlap and the submodular function will be further from being
#monotone.

N=100;
U_incid = Array{Array}(N);
for i=1:N
	U_incid[i] = []
end
for i = 1:N
		degree_of_i = 0.15*N + Int(round(rand()*N/2));
		for j = 1:degree_of_i
			random_neighbor = Int(round(rand()*(N-1)+1));
			push!(U_incid[i], random_neighbor);
		end
end

#assign random weights to vertices on both sides of the bipartition
m = 100*rand(N,1); #for vertices in U
weight = 100*rand(N, 1); #for vertices in V

### here is our very terrible submodular function
function jegelka_example(S)
	val = sum(m[S]);
	neighborhood = [];
	for i in S
			neighborhood = vcat(neighborhood, U_incid[i]) #compute the neighborhood of the selection of vertices from U
	end
	neigborhood = unique(neighborhood);
	val = val + 10*sqrt(sum(weight[neighborhood])); # one can vary the constant in front of the weight of the neighborhood as well.
	return val;
end

### now construct an optimization problem over the associated submodular polytope
a = rand(N)
### ??? I have no idea what the right syntax here is
g = DiffableFunction(x -> 1/2*sum((x-a).^2))
constraint(S) = AssocPoly(jegelka_example(S))
problem = CombiProblem(g, constraint)
### end ???

x = frank_wolfe_away(problem, verbose=true, maxiters=10)
