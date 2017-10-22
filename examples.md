# Best min-norm-point algos in practice

* Fugishige-Wolfe
* FW with away steps
* Projected subgradient

No implementable strongly polynomial algorithm is < O( SFO * n^5 + n^6)

# Convex functions

We want examples that are smooth and nonsmooth, strongly and non-strongly convex,
and polyhedral or SDP-based.

1. Least squares (strongly convex and smooth)

2. Huber regression (not strongly convex but smooth)

3. L1 regression / total variation (neither strongly convex nor smooth)

4. Logistic regression (and maybe Multinomial regression)

5. Hinge loss classification

6. Hinge loss classification + l2 regularization (strongly convex but not smooth)

7. Log determinant

8. Lambdamax

# Submodular functions

1. Cardinality based functions

2. Weighted graph cuts

3. Jegelka crazy example

4. Log determinant

  * Determinental point processes - subsets sampled ~ volume of subset
  * Submodular function selects indices
  * See Djolonga and Krause
