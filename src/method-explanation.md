This formula is more sophisticated than the procedure you originally described. It is not evaluating individual variables directly; it is evaluating a *subset of variables* according to how well that subset reproduces the PCA structure of the full dataset.

The criterion is:

\[ r_m = \sqrt{
\frac{\sum_{i=1}^{p}\lambda_i r_{mi}^2}
{\sum_{j=1}^{p}\lambda_j}
}\]

where:

- (\lambda\_i) is the eigenvalue of PC (i) from the full dataset;
- (r\_{mi}\^2) is the squared multiple correlation between PC (i) and the subset of variables;
- 
  (m) indexes a candidate subset.

The numerator computes a weighted average of how well the subset predicts each principal component, with the weights given by the variance explained by that component ((\lambda\_i)). The denominator normalizes by total variance.

Conceptually, this means:

> "How much of the overall multivariate information contained in the full PCA can be recovered from this subset of variables?"

That is a very sensible objective if the goal is dimension reduction.

A few observations:

1.  The method is fundamentally a subset-selection method

It does not identify variables that are individually important.

Instead it asks:

> Among all subsets of size (k), which one best reproduces the PCA solution?

This is closer to "best subset selection" than to a variable importance measure.

2.  It has a clear interpretation

The paper's statement that (r_m\^2) can be interpreted as the percentage of total variation explained by the subset is essentially correct.

If:

\[

r_m\^2 = 0.90

\]

then the selected variables capture about 90% of the information represented by the full PCA solution.

This is attractive because the metric is interpretable.

3.  It avoids one problem of your original proposal

Earlier I noted that simply weighting loadings by eigenvalues can favour redundant variables.

Suppose variables A, B, and C are almost identical.

A variable-ranking approach might select all three because they load strongly on PC1.

This subset criterion tends to be less vulnerable because once A is included, adding B and C contributes relatively little additional ability to reconstruct the PCs.

So it has some built-in protection against redundancy.

4.  It assumes PCA defines "information"

The largest conceptual limitation remains.

The method equates information with the variance structure captured by PCA.

This is appropriate if your objective is statistical data reduction.

It may be less appropriate if your index is intended to represent a substantive concept.

For example, in a poverty index, a variable with low variance but high policy relevance could be discarded because it contributes little to reproducing the PCA.

5.  There is a connection to PCA biplots and ordination preservation

The paper mentions matching the "principal component ordination."

This is important.

The criterion is essentially asking whether observations would occupy approximately the same positions in PCA space if only the subset of variables were used.

In ecology and multivariate statistics this is often a desirable property.

6.  Computational issues

The challenge is that finding the optimal subset of size (k) requires searching over

\[\binom{p}{k}\]

possible subsets.

For 30 variables and a subset size of 10, that is over 30 million combinations.

Most implementations therefore use heuristic searches (stepwise procedures, branch-and-bound algorithms, genetic algorithms, etc.).

For index construction, I think this criterion is substantially stronger than your original variable-level weighting idea. Your original idea measures how strongly a variable aligns with dominant PCs. The paper's criterion measures how well a set of variables reproduces the entire multivariate structure of the data. Statistically, that is a much more coherent objective for indicator reduction.

The key question is whether your index aims to preserve the data's covariance structure (where this method is attractive) or to represent a theoretical construct (where substantive considerations should still constrain variable selection).