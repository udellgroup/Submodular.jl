using HDF5, JLD
using Distributions

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

m = 100*rand(N,1);
weight = 100*rand(N, N);


function jegelka_example(S)
	val = 0;
	for i in S 
		val = val + m(i);
		w = 0;
		for j in U_incid(i)
			w = w + weight(i, j);
		end
		val = val + sqrt(w);
	end
	return val;
end

