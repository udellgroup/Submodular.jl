#Function to minimize Bregman divergences h() with respect to a given 
#vector y over the submodular base polytope of a cardinality-based submodular function 

using Distributions

function cardincfix(g, y_real, divergence) 
  # here F(S) = g(|S|) for some non-decreasing concave function g, g(0)=0.
  # currently implemented for divergences: Euclidean
  # check if g is concave, non-decreasing, g is an Array for now
  
    N = length(y_real);
	if(divergence=="euclidean")
		f_obj(x,y) = norm(x-y)^2;
	end
	
	function f_submod(k)
		if(k==0) return 0 
		else 
			return g[k]	
		end
	end
			
	sorted_indices = sortperm(y_real, rev=true);
	y = y_real[sorted_indices]
	if(y[1]<=0 && divergence!="euclidean")
   		error("y out of range for itakura-saito, logistic, entropy");
	elseif(y[1]>=1&&divergence=="logistic")
    	error("y out of range for logistic");
	end
	
	counter = 0; #maintains the tight set
	epsilon = zeros(N,1); 

	if(divergence=="euclidean")
    	x_opt = zeros(N,1);
    	x_opt = round.(x_opt, 4);
    	feas_x_opt = (g[N]/N)*ones(N,1);
    	while(counter<N) 
        	for j = counter+1:N
            	if(counter==0)
                	epsilon[j] = (f_submod(j) - sum(y[1:j]))/j;
            	else
                	epsilon[j] = (f_submod(j) - f_submod(counter) - sum(y[counter+1:j]))/(j-counter);
            	end
        	end
        	(e_k, I) = findmin(epsilon[counter+1:N]);
        	k = I + counter;   
        	x_opt[counter+1:k] = y[counter+1: k] + e_k;
        	counter = k;
        	temp = counter+1;
        	while(temp<=N && x_opt[temp] - y[temp] <= e_k)
            	x_opt[temp] = e_k+ y[temp];
            	temp = temp+1;
        	end
    	end
	end
	temp = zeros(N,1);
	for i =1:N
		temp[sorted_indices[i]] = x_opt[i]; 
	end
	return temp;
end


################################
#example: 
N=100;
g = rand(N);
g = sort(g, rev=true);
g = (g + 0.1)/1.01;
for i = 2: N
	g[i] = g[i] + g[i-1];
end

# g is a concave non-decreasing function

y = rand(N); # point to project on the base polytope of F(S) = g(|S|)

euclidean_proj = cardincfix(g, y, "euclidean");

#sanity checks: 
sortperm(y) == sortperm(euclidean_proj[1:N])
sum(euclidean_proj)>=g[N]-0.0001 ## sometimes there is an error of 10^-10 in these two numbers. 
sum(euclidean_proj)<=g[N]+0.0001

##############################

