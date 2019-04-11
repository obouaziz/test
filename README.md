# test

Implementation of corrected two sample tests: corrections for the Pearson, Kendall and Spearman correlation tests, 
 the Mann-Whitney (Wilcoxon) rank sum test, the Mann-Whitney (Wilcoxon) signed rank test and a variance test.
 The package also proposes a test for the median based on a kernel estimator of the density and a confidence interval for the median based on rank statistics. All the tests are asymptotically calibrated meaning that
the probability of rejection under the null hypothesis is asymptotically equal to 5%.

The package only contains the functions:

- `cortest` that implements the corrected Pearson, Kendall's and Spearman's tests.
- `wilcoxtest` that implements the corrected version of the Wilcoxon (Mann-Whitney) test for two independent samples and for two paired samples.
- `vartest` that implements a variance test based on the Welch correction for the variables (x_i-mean(x))^2 and (y_i-mean(y))^2.
- `mediantest` that implements a median test based on a kernel estimator of the density at 0.
- `medianCI` that produces a confidence interval for the median based on the rank statistic.
- `tiebreak` which randomly breaks ties in vectors, either inside the vector or between two vectors.
