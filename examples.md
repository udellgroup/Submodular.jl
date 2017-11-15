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

## Sample code to generate loss functions
```
    n =
    m =  # probably choose m = 2n for now; choose m < n for non-unique solution
    x♮ = randn(n)
    A = randn(m,n)

    ### regression losses
    b = A*x♮ + .1*randn(m)

    # least squares loss
    fquad(x) = norm(A*x - b)^2

    # huber loss
    fhuber(x) = sum(huber.(A*x - b))

    # l1 loss
    fquad(x) = norm(A*x - b, 1)

    ### classification losses
    boolb = sign(A*x♮ + .1*randn(m))

    # logistic loss
    f(x) = sum(log.(1+exp.(-boolb.*A*x)))

    # hinge loss
    f(x) = sum(max.(1-boolb.*A*x, 0))

    # hinge loss + l2
    f(x) = sum(max.(1-boolb.*A*x, 0)) + norm(x)^2
```

# Submodular functions

1. Cardinality based functions

2. Weighted graph cuts

3. Jegelka crazy example

4. Log determinant

  * Determinental point processes - subsets sampled ~ volume of subset
  * Submodular function selects indices
  * See Djolonga and Krause
