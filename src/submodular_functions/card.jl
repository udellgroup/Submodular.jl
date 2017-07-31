# Handles functions of the form f(A) = g(length(A)), where g is a concave function

type Card <: Function
  child::Function
end

function in(w::AbstractVector, p::SubPoly)
  n = length(w)
  @assert length(V) == n
  checker = true
  S = sort(V)
  ordering = sortperm(w, rev = true)
  h = 0
  ind = []
  for i = 1:length(w)
    h += - f(S[ind]) - w[i] + f(S[push!(ind, ordering[i])])
    if h < 0
      checker = false
      break
    end
  end
  return checker
en
