# test

Implementation of corrected two sample tests: corrections for the Pearson, Kendall and Spearman correlation tests, 
 the Mann-Whitney (Wilcoxon) rank sum test, the Mann-Whitney (Wilcoxon) signed rank test and a variance test. 
 The package also proposes a test for the median based on a kernel estimator of the density and a confidence interval for the median based on rank statistics. All the tests are asymptotically calibrated meaning that
the probability of rejection under the null hypothesis is asymptotically equal to 5%. The package also proposes a test for independence between two continuous variables, based on the maximum distance between the joint empirical cumulative distribution function and the product of the marginal
empirical cumulative distribution functions is also implemented. The distribution of this test has been numerically obtained and this is also an exact test.

The package contains the functions:

- `cortest` that implements the corrected Pearson, Kendall's and Spearman's tests. As compared to the original tests in cor.test which all assume independence between the variables under the null hypothesis, the corrected tests assume that the correlation (of the different types) is equal to 0 under the null.
- `wilcoxtest` that implements the corrected version of the Wilcoxon (Mann-Whitney) test for two independent samples and for two paired samples. In the two independent sample case, as compared to the original test in wilcox.test which assumes that the variables have the same distribution under the null hypothesis, the corrected test assumes under the null that the probability that one variable exceeds the other is equal to 0.5.
- `vartest` that implements a variance test based on the Welch correction for the variables (x_i-mean(x))^2 and (y_i-mean(y))^2. As compared to the original test in var.test which only works under the gaussian scenario, the corrected test works for any distribution of the two variables as long as the fourth order moments exist for both variables.
- `mediantest` that tests if the median of the random variable (one sample case) or of the difference between two random variables (two sample case) is equal to 0.  The test is based on the asymptotic law of the empirical median and uses a kernel estimator to estimate the density of the variable (in the one sample case) or of the difference between the two variables in the two sample case.
- `medianCI` that produces a confidence interval for the median of the random variable (one sample case) or of the difference between two random variables (two sample case) based on the rank statistic of the median.
- `robustest`that tests the independence between two continuous variables. The test is based on the maximum distance between the joint empirical cumulative distribution function and the product of the marginals. The distribution of this test has been numerically obtained, the test is exact for all n<=150 and approximated for n>150.
- `tiebreak` which randomly breaks ties in vectors, either inside the vector or between two vectors.

The dataset `Evans` can also be loaded from the **test** package. A description of this dataset can be found in the **lbreg** package.
