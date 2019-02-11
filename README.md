# test

Implementation of corrected two sample tests. A corrected version of the Pearson and Kendall correlation tests, 
 the Mann-Whitney (Wilcoxon) rank sum test, the Mann-Whitney (Wilcoxon) signed rank test and the Fisher variance test are implemented.
 The package also proposes a test for the median. All these corrected tests are asymptotically calibrated meaning that
the probability of rejection under the null hypothesis is asymptotically equal to 5%.

For the moment the package only contains the function:

- \code{cortest} that implements the corrected Pearson and Kendall's tests.
- \code{wilcoxtest} that implements the corrected version of the Wilcoxon (Mann-Whitney) test for two independent samples and for two paired samples.
