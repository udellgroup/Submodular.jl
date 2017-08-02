# computes the proximal operator of the indicator functino of the base polyhedron
function prox(p::BasePoly, w::AbstractVector)
  # step 1 initialization
  n = length(w)
  V = sort(p.V)     # the sorted base set
  x = zeros(n)      # the candidate
  for i = 1:n
    x[i] = f(V[1:i]) - f(V[1:i-1])
  end
  S = x             # the set of extreme points that constitute the affine set
  λ = [1]           # the coefficients of extreme points

  # step 2 find the point z that minimize (x-w)^(z-w) st z in the base polyhedron
  z = fenchel(p, x - w)
  minprod = dot(z - w, x - w)
  curprod = dot(x - w, x - w)

end
