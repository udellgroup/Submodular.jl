#Function to minimize Bregman divergences h() with respect to a given
#vector y over the submodular base polytope of a cardinality-based submodular function

using Distributions

export card_inc_fix

function card_inc_fix(F::CardBasedAtom,
	                    w₀::AbstractArray,
										  divergence)
  # here F(S) = g(|S|) for some non-decreasing concave function g, g(0)=0.
  # currently implemented for divergences: Euclidean
  # check if g is concave, non-decreasing, g is an Array for now

  n = length(w₀)
  g = zeros(n)
  for i = 1:n
    g[i] = evaluate(F, collect(1:i))
  end

	# if(divergence=="euclidean")
	# 	f_obj(x,y) = norm(x-y)^2
	# end
	#
	# function f_submod(k)
	# 	if(k==0) return 0
	# 	else
	# 		return g[k]  # the submodular function
	# 	end
	# end

	sorted_indices = sortperm(w₀, rev=true)
	w = w₀[sorted_indices]
	if(w[1]<=0 && divergence!="euclidean")
   	error("The input is out of range for itakura-saito, logistic, entropy")
	elseif(w[1]>=1&&divergence=="logistic")
  	error("The input is out of range for logistic")
	end

	counter = 0                             # maintains the tight set
	epsilon = zeros(n,1)

	if(divergence == "euclidean")
    x_opt = zeros(n,1)
  	# x_opt = round.(x_opt, 4)
  	feas_x_opt = (g[n]/n) * ones(n,1)
  	while(counter < n)
      for j = (counter + 1) : n
        if(counter == 0)
          epsilon[j] = (g[j] - sum(w[1:j]))/j
        else
          epsilon[j] = (g[j] - g[counter] - sum(w[counter+1:j]))/(j-counter)
        end
      end
      (e_k, I) = findmin(epsilon[counter+1:n])
      k = I + counter
    	x_opt[counter+1:k] = w[counter+1: k] + e_k
      counter = k
    	temp = counter+1
      while(temp<=n && x_opt[temp] - w[temp] <= e_k)
        x_opt[temp] = e_k+ w[temp]
        temp = temp+1
      end
    end
	end
	temp = zeros(n,1)
	for i =1:n
		temp[sorted_indices[i]] = x_opt[i]
	end
	return temp
end

################################

n = 20
println("Example on a ground set of $n elements:")

g = rand(n)
g = sort(g, rev=true)
g = (g + 0.1)/1.01
h = copy(g)
for i = 2: n
  g[i] = g[i] + g[i-1]
end
# g is a concave non-decreasing function
println("Cardinality-based submodular function F(S) = g(|S|), where g is given by:")
println("g = ",g)
y = rand(n) # point to project on the base polytope of F(S) = g(|S|)
println()
println("The random point y that we are projecting on the base polytope B(F) is:")
println("y = ", y)

S = SetVariable(n)
gg = card(h, S)


euclidean_proj = card_inc_fix(gg, y, "euclidean")
println()
println("The computed Euclidean projection is: ")
println("projection = ", euclidean_proj)

#sanity checks:
println()
println("Performing sanity checks: ")
println("Is the sorted order of indices in y and the projection the same?")
#this is a known property of projections under uniform divergences over cardinality-based polytopes
println(sortperm(y) == sortperm(euclidean_proj[1:n]))

println("Is the constraint x(E) = F(E) = g(n) satisfied upto an error of 10^-10?")
sum(euclidean_proj)>=g[n]-(10.0)^(-10)
sum(euclidean_proj)<=g[n]+(10.0)^(-10)

##############################
