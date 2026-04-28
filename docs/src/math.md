# Survey Sampling Math

Below, you'll find informal derivations of the formulas used in the library. 

## Sampling Without Replacement

- Let the indicator variable ``I_i`` be ``1`` if individual ``i`` is
  included in the sample. Let ``\pi_i = E[I_i]``.

- Let ``\breve{y}_i = y_i / \pi_i``. In vector form, ``\breve{y}`` is a
  size ``n`` vector of the ``\breve{y}_i`` in the sample.

- An unbiased sample sum estimator is
  ``\sum_{i \in U} I_i \breve{y}_i``.

- Variance of this estimator is
  ``\sum_{i,j \in U} \Delta_{ij}\breve{y}_i \breve{y}_j``. where
  ``\Delta_{i,j} = \text{Cov}(I_i, I_j)``

- An estimator of the variance is therefore
  ``\sum_{i,j \in U} I_i I_j \breve{\Delta}_{ij} \breve{y}_i\breve{y}_j``
  where ``\breve{\Delta}_{ij} = \frac{\Delta_{ij}}{E[I_iI_j]}``.

- In matrix form, that's ``\breve{y} \breve{\Delta} \breve{y}``.

## With Replacement

- Let ``p_i`` be the probability of sampling individual ``i`` at each
  round.
- Sum estimator is ``\frac{1}{n}\sum_i \frac{y_i}{p_i}``.
- Variance estimator is ``\text{Var}(\frac{y_i}{p_i})``.

## General Taylor Series

- Say we don't observe ``y_i = f(z_i)`` but we do observe ``z_i``. Let
  ``a = \sum_{i \in U} z_i``.
- Sum estimator is ``f(\hat{a})`` where ``\hat{a}`` is the vector of
  component sum estimates.
- ``f(\hat{a}) \approx f(a) + \nabla f(a)(\hat{a} - a)``.
- So
  ``\text{Var}(f(\hat{a})) \approx \text{Var}(\langle \nabla f(a), \hat{a}\rangle)``
- We can use ``\nabla f(\hat{a})`` to estimate ``\nabla f(a)``.

## Taylor Series Variance Without Replacement

Plug in sum estimator to get
``\text{Var}(\nabla f(a)^T \sum_i \breve{z}_i)``. This is the variance
of a different sum, giving ``g^T \breve{\Delta}g`` where
``g = \langle \nabla f(T), \breve{z} \rangle``.

## Taylor Series Variance With Replacement

Plug in sum estimator to get
``\text{Var}(\nabla f(a)^T \frac{1}{n}\sum_i \breve{z}_i)``. This is the
variance of a different sum, giving ``\text{Var}(g)`` where
``g = \langle \nabla f(T), \breve{z} \rangle``.

## Regression Estimates

- Goal: population regression coefficient ``\beta = T^{-1}t`` where
  ``T = \sum_{i \in U} x_i x_i^T`` and ``t=\sum_{i \in U} x_i y_i``.
- We know ``s = \sum_{i \in U} x_i``. So we can estimate the total with
  ``f(T, t) = s^T\beta`` using sampling estimates ``\hat{T}`` and
  ``\hat{t}``.
- For variance, use Taylor approximation and drop the constant
  intercept:
  ``\text{Var}(f(\hat{T}, \hat{t})) \approx \text{Var}( \langle \nabla f(T, t), \begin{bmatrix} \hat{T} - T & \hat{t} - t \end{bmatrix} \rangle )``
- Matrix calculus gives ``\nabla_t f = s^T T^{-1}`` and
  ``\nabla_T f= -T^{-1}(st^T)T^{-1}``

```math

\begin{align*}
\langle \nabla f(T, t), \begin{bmatrix} \hat{T} - T & \hat{t} - t \end{bmatrix} \rangle &= s^T T^{-1}(\hat{t} - t) - s^T T^{-1} (\hat{T} - T)T^{-1}t\\
&= s^T T^{-1} ((\hat{t} - t) - (\hat{T} - T)\beta )\\
&= s^T T^{-1} ((\hat{t} - T\beta) - (\hat{T} - T)\beta)\\
&= s^T T^{-1}(\hat{t} - \hat{T}\beta)\\
\end{align*}

```

## Regression Without Replacement

- We can use the sampling approximations ``\hat{T} = X^T\Pi^{-1}X`` to
  ``T = \sum_{i \in U} x_i x_i^T`` and ``\hat{t} =X^T\Pi^{-1}y`` to
  ``t=\sum_{i \in U} x_i y_i``.

- Substitute this into the variance expression to get
  ``\text{Var}(f) = s^T T^{-1}X^T \Pi^{-1}(y - X \beta)``

- This means we're calculating the variance of a without-replacement
  population sum: ``\frac{g_i}{\pi_i}`` where
  ``g_i = s^T T^{-1} x_i(y_i - x_i^T\beta)``. We can use ``\hat{T}`` to
  approximate ``T`` here.

- We know that's ``g^T \breve{\Delta} g``.

## Regression With Replacement

- We can use the sampling approximations ``\hat{T} = X^T\Pi^{-1}X`` to
  ``T = \sum_{i \in U} x_i x_i^T`` and ``\hat{t} =X^T\Pi^{-1}y`` to
  ``t=\sum_{i \in U} x_i y_i`` where ``\Pi_{ii} = np_i``
- This means we're calculating the variance of a with-replacement
  population sum: ``\frac{g_i}{\pi_i}`` where
  ``g_i = T^{-1} x_i(y_i - x_i^T\beta)``
- That's ``\text{Var}(g_i)``.

## Sampling Clusters Without Replacement

We can apply the without-replacement estimator recursively to estimates of cluster totals.
In this context, let $I_i = 1$ if cluster $i$ is sampled. By linearity, the estimator ``\hat{t} = \sum_{i \in \text{PSU}} I_i \frac{\hat{t}_i}{\pi_i}`` is unbiased. We can get its variance from the law of total variance:

```math
\begin{align*}
\text{Var}(\hat{t}) &= \text{Var}(E[\hat{t} | I]) + E[\text{Var}(\hat{t} | I)] \\
&= V_I\left[\sum_{i \in \text{PSU}} I_i \frac{\hat{t}}{\pi_i}\right] + E_I\left[\sum_{i \in \text{PSU}} I_i \frac{\text{Var}(\hat{t})}{\pi_i}\right]
\end{align*}
```

The first term is the variance of the without-replacement estimator over ``\hat{t}_i`` values. The second term is the with-replacement estimator of the sum of the cluster variances. We already know how to compute each of these!

## Taylor Series Estimates with Clustering

When we're estimating a function of population totals ``\hat{t} = f(\hat{a})`` where the population totals are aggregated from cluster totals, we can estimate the variance in a similar fashion.

```math
\begin{align*}
\text{Var}(E[\hat{t} | I]) &= \text{Var}_I(\langle \nabla f, \hat{a} \rangle) \\
&\approx \text{Var}_I(\langle \nabla f, \sum_{i \in \text{PSU}} I_i \frac{\hat{z}_i}{\pi_i} \rangle) \\
&= \text{Var}_I(\sum_{i \in \text{PSU}}I_i \langle \nabla f, \frac{\hat{z}_i}{\pi_i} \rangle)
\end{align*}
```

This is the variance of the without-replacement estimator over ``\langle \nabla f, \hat{z}_i`` values.

```math
\begin{align*}
E[\text{Var}(\hat{t} | I)] &= E[\text{Var}(f(\hat{a} | I)] \\ 
&\approx E[\text{Var}(\langle \nabla f, \hat{a} \rangle | I)] \\ 
&= E\left[\sum_{i \in \text{PSU}} \frac{I_i}{\pi_i} (\langle \nabla f, \text{Var}(\hat{z}_i) \rangle)\right] \\ 
\end{align*}
```
This is the without-replacement estimator of the sum of ``\langle \nabla f, \text{Var}(\hat{z}_i)`` values.
