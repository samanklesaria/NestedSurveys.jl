# Survey Sampling Math

Below, you'll find informal derivations of the formulas used in the
library.

## The Design-Based Perspective

Most statistical inference assumes something about the data-generating
process. Linear regression assumes errors are normally distributed.
Bayesian methods assume a prior. But survey statisticians have a
different trick available: *they control the sampling process*, so they
know its probability distribution exactly.

This is the **design-based** (or **randomization-based**) approach.
Instead of assuming the population follows some model, you treat the
population values as fixed unknowns and derive your inference from the
randomness in *how you sampled*. Your uncertainty isn't about the world
--- it's about which people you happened to pick.

The key object is the **inclusion probability** ``\pi_i``: the
probability that individual ``i`` appears in your sample. If you know
these, you can do a lot.

## The Horvitz-Thompson Estimator

Let the indicator variable ``I_i`` be ``1`` if individual ``i`` is
included in the sample, so ``\pi_i = E[I_i]``. Let
``\breve{y}_i = y_i / \pi_i``. In vector form, ``\breve{y}`` is a size
``n`` vector of the ``\breve{y}_i`` in the sample.

An unbiased sample sum estimator is
$$\hat{t} = \sum_{i \in U} I_i \breve{y}_i$$
This is known as the *Horvitz-Thompson estimator*. The intuition: if
someone had a 10% chance of being sampled and showed up, they probably
"represent" about 10 people in the population, so we weight their value
by ``1/\pi_i = 10``.

The variance is where things get slightly more involved. Define
``\Delta_{ij} = \text{Cov}(I_i, I_j)``. Then:

$$\text{Var}(\hat{t}) = \sum_{i,j \in U} \Delta_{ij} \breve{y}_i \breve{y}_j = \breve{y}^T \Delta \breve{y}$$
To estimate this variance from the sample, we divide each
``\Delta_{ij}`` by the joint inclusion probability ``E[I_i I_j]``,
giving ``\breve{\Delta}_{ij}``. The estimated variance is then
``\breve{y}^T \breve{\Delta} \breve{y}`` --- a sum only over observed
units.

### Simple Random Sampling

For simple random sampling without replacement (SRS) from a population
of ``N``, drawing ``n`` units:

- ``\pi_i = n/N`` for all ``i``
- ``\pi_{ij} = \frac{n(n-1)}{N(N-1)}`` for ``i \neq j``
- ``\Delta_{ij} = \pi_{ij} - \pi_i \pi_j``

The variance simplifies to the familiar ``N^2 \cdot \frac{1-f}{n} S^2``
where ``f = n/N`` is the sampling fraction and ``S^2`` is the population
variance --- adjusted by the **finite population correction** ``(1-f)``.
When you sample a large fraction of the population, your estimate gets
better, and this correction captures that.

## Sampling With Replacement

When units are drawn independently with probability ``p_i`` at each
draw, the *Hansen-Hurwitz* estimator is:

$$\hat{t} = \frac{1}{n} \sum_{i=1}^{n} \frac{y_i}{p_i}$$
Each draw yields an unbiased estimate of the total, and you average
``n`` of them. The variance estimator is just the sample variance of the
``y_i / p_i`` values. This is simpler to compute because the draws are
independent --- you don't need joint inclusion probabilities.

## Beyond Totals: Taylor Series Linearization

Most quantities of interest aren't linear in the population values.
Means, ratios, regression coefficients --- they're all nonlinear
functions of sums. The standard tool for handling this is **Taylor
series linearization**. Say we observe want to compute ``f(a)`` where
``a = \sum_{i \in U} z_i`` is a vector of population totals of different
variables and ``z_i`` are the observation vectors from our survey. As
long as ``f`` is continuous, a consistent estimator is ``f(\hat{a})``
where ``\hat{a}`` is the vector of component sum estimates.

To get the variance of this estimator, note that
``f(\hat{a}) \approx f(a) + \nabla f(a)(\hat{a} - a)``. So
``\text{Var}(f(\hat{a})) \approx \text{Var}(\langle \nabla f(a), \hat{a}\rangle)``.
We can use ``\nabla f(\hat{a})`` to estimate ``\nabla f(a)``.

If the underlying ``z_i`` were sampled without replacement, we could
plug in the Horvitz-Thompson estimator as ``\hat{a}`` to get

```math

\text{Var}(f(\hat{a})) \approx \text{Var}\left(\nabla f(a)^T \sum_i \breve{z}_i\right)

```

This is the variance of a different sum, giving ``g^T \breve{\Delta}g``
where ``g = \langle \nabla f(T), \breve{z} \rangle``.

If the ``z_i`` were sampled *with* replacement, we could plug in the
Hansen-Hurwitz estimator to get

$$\text{Var}(f(\hat{a})) = \text{Var}\left(\nabla f(a)^T \frac{1}{n}\sum_i \frac{y_i}{p_i}\right)$$
Once again, this is the variance of a sum, which we already know how to
compute!

## Regression-Assisted Estimation

Sometimes you have auxiliary information: population totals of
covariates that aren't available at the unit level in the sample.
Regression-assisted estimation exploits this.

You want to estimate ``s^T \beta`` where ``\beta = T^{-1} t`` is the
population regression coefficient, ``T = \sum_{i \in U} x_i x_i^T``,
``t = \sum_{i \in U} x_i y_i``, and ``s = \sum_{i \in U} x_i`` is the
known population covariate total.

- We know ``s = \sum_{i \in U} x_i``. So we can estimate the total with
  ``f(T, t) = s^T\beta`` using sampling estimates ``\hat{T}`` and
  ``\hat{t}``.
- For variance, use Taylor approximation and drop the constant
  intercept:

```math

\text{Var}(f(\hat{T}, \hat{t})) \approx \text{Var}( \langle \nabla f(T, t), \begin{bmatrix} \hat{T} - T & \hat{t} - t \end{bmatrix} \rangle )

```

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

Simplifying further requires us to think about which estimators we're
using for ``\hat{T}`` and ``\hat{t}``.

## Regression Without Replacement

If we've been sampling *without* replacement, we can use the sampling
approximations ``\hat{T} = X^T\Pi^{-1}X`` to
``T = \sum_{i \in U} x_i x_i^T`` and ``\hat{t} =X^T\Pi^{-1}y`` to
``t=\sum_{i \in U} x_i y_i``. Substitute this into the variance
expression to get

```math

\text{Var}(f) = \text{Var}(s^T T^{-1}X^T \Pi^{-1}(y - X \beta))

```

The term in the ``\text{Var}`` expression is just a Horvitz-Thompson sum
estimate of population values
``g_i = s^T T^{-1} x_i(y_i - x_i^T\beta)``. We know that's
``g^T \breve{\Delta} g``. We can use ``\hat{T}`` to approximate ``T``
here.

## Regression With Replacement

If we've been sampling *with* replacement, we can use the sampling
approximations ``\hat{T} = X^TP^{-1}X`` to
``T = \sum_{i \in U} x_i x_i^T`` and ``\hat{t} =X^TP^{-1}y`` to
``t=\sum_{i \in U} x_i y_i`` where ``P_{ii} = np_i``.

If you squint a bit, you'll see that this is just the variance of the
Hansen-Hurwitz sum estimator on population values
``h_i = s^TT^{-1} x_i(y_i - x_i^T\beta)``.

## Multi-Stage Cluster Sampling

Here's where things get genuinely interesting, and where most survey
software starts to strain.

Real surveys rarely sample individuals directly. You might first sample
counties (making them *primary sampling units* or PSUs), then sample
households within those counties. This is two-stage cluster sampling,
and the nesting can go arbitrarily deep: counties → blocks → households
→ people.

Let ``I_i = 1`` if cluster ``i`` is selected at the first stage, and let
``\hat{t}_i`` be the estimated total within cluster ``i`` from the
second stage. The overall estimator is:

$$\hat{t} = \sum_{i \in \text{PSU}} I_i \frac{\hat{t}_i}{\pi_i}$$
This is unbiased by linearity. For its variance, apply the law of total
variance:

```math

\begin{align*}
\text{Var}(\hat{t}) &= \text{Var}(E[\hat{t} | I]) + E[\text{Var}(\hat{t} | I)] \\
&= V_I\left[\sum_{i \in \text{PSU}} I_i \frac{\hat{t}}{\pi_i}\right] + E_I\left[\sum_{i \in \text{PSU}} I_i \frac{\text{Var}(\hat{t})}{\pi_i}\right]
\end{align*}

```

The first term is the variance of the without-replacement estimator over
``\hat{t}_i`` values. The second term is the with-replacement estimator
of the sum of the cluster variances. We already know how to compute each
of these!

This recursive structure extends naturally to any depth of nesting.
Three-stage sampling? The within-cluster variance is itself computed via
a two-stage formula.

## Taylor Series Estimates with Clustering

For nonlinear estimators like ratios, the same decomposition applies via
Taylor linearization. The gradient ``\nabla f`` filters through to
weight the cluster-level sums, and everything falls out cleanly.

```math

\begin{align*}
\text{Var}(E[\hat{t} | I]) &= \text{Var}_I(\langle \nabla f, \hat{a} \rangle) \\
&\approx \text{Var}_I(\langle \nabla f, \sum_{i \in \text{PSU}} I_i \frac{\hat{z}_i}{\pi_i} \rangle) \\
&= \text{Var}_I(\sum_{i \in \text{PSU}}I_i \langle \nabla f, \frac{\hat{z}_i}{\pi_i} \rangle)
\end{align*}

```

This is the variance of the without-replacement estimator over
``\langle \nabla f, \hat{z}_i`` values.

```math

\begin{align*}
E[\text{Var}(\hat{t} | I)] &= E[\text{Var}(f(\hat{a} | I)] \\ 
&\approx E[\text{Var}(\langle \nabla f, \hat{a} \rangle | I)] \\ 
&= E\left[\sum_{i \in \text{PSU}} \frac{I_i}{\pi_i} (\langle \nabla f, \text{Var}(\hat{z}_i) \rangle)\right] \\ 
\end{align*}

```

This is the without-replacement estimator of the sum of
``\langle \nabla f, \text{Var}(\hat{z}_i)`` values.
