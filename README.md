# test

Implementation of corrected two sample tests: corrections for the Pearson, Kendall and Spearman correlation tests, 
 the Mann-Whitney (Wilcoxon) rank sum test, the Mann-Whitney (Wilcoxon) signed rank test and a variance test. 
 The package also proposes a test for the median based on a kernel estimator of the density and a confidence interval for the median based on rank statistics. All the tests are asymptotically calibrated meaning that
the probability of rejection under the null hypothesis is asymptotically equal to 5%. The package also proposes a test for independence between two continuous variables, based on the maximum distance between the joint empirical cumulative distribution function and the product of the marginal
empirical cumulative distribution functions is also implemented. The distribution of this test has been numerically obtained and this is also an exact test.

The package contains the functions:

- `cortest` that implements the corrected Pearson, Kendall's and Spearman's tests. As compared to the original tests in cor.test which all assume independence between the variables under the null hypothesis, the corrected tests assume that the correlation (of the different types) is equal to 0 under the null.
- `wilcoxtest` that implements the corrected version of the Wilcoxon (Mann-Whitney) test for two independent samples and for two paired samples. In the two independent samples case, as compared to the original test in wilcox.test which assumes that the variables have the same distribution under the null hypothesis, the corrected test assumes that the probability that one variable exceeds the other is equal to 0.5 under the null.
- `vartest` that implements a variance test based on the Welch correction for the variables (x_i-mean(x))^2 and (y_i-mean(y))^2. As compared to the original test in var.test which only works under the gaussian scenario, the corrected test works for any distribution for the two variables as long as the fourth order moments exist for both variables.
- `medianCI` that produces a confidence interval for the median based on the rank statistic.
- `tiebreak` which randomly breaks ties in vectors, either inside the vector or between two vectors.

The dataset `Evans` can also be loaded from the **test** package. A description of this dataset can be found in the **lbreg** package.
