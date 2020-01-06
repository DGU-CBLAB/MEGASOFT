## METASOFT C++

## Table of Content
- [About](#about)
- [How to Use](#use)
- [Prerequisite](#prerequisite)
- [Reference](#reference)

#
## About <a name="about"></a>
In Genome-wide association Study Analysis(GWAS), <a href="http://genetics.cs.ucla.edu/meta/">METASOFT</a> is considered a widely used software, which is effective in performing wide range of basic and complicated meta-analytic methods. However, <a href="http://genetics.cs.ucla.edu/meta/">METASOFT</a> can be improved in terms of computing performance. This repository is another version of <a href="http://genetics.cs.ucla.edu/meta/">METASOFT</a> where computing load balancing is used to improve computing time.

#

## How to Use<a name="use"></a>

<a href="http://genetics.cs.ucla.edu/meta/">Overview Link</a>

User's Guide (differ from original METASOFT)
```
Metasoft.exe [option]
-input <FILE>                       Input file (Required)
    -output <FILE>                      Output file (default='out')
    -log <FILE>                         Log file (default='log')
    -pvalue_table <FILE>                Pvalue table file (default='HanEskinPvalueTable.txt')
    -lambda_mean <FLOAT>                (Random Effects) User-specified lambda for mean effect part
                                        (default=1.0)
    -lambda_hetero <FLOAT>              (Random Effects) User-specified lambda for heterogeneity
                                        part (default=1.0)
    -binary_effects                     Compute binary effects model p-value (default=false)
    -binary_effects_sample <INT>        (Binary effects) Number of importance sampling samples
                                        (default=1,000)
    -binary_effects_p_thres <FLOAT>     (Binary effects) P-value threshold determining if we will
                                        use large number of samples (default=1E-4)
    -binary_effects_large <INT>         (Binary effects) Large number of importance sampling samples
                                        for p-values above threshold (default=100,000)
    -mvalue                             Compute m-value (default=false)
    -mvalue_method <METHOD_NAME>        Which method to use to calculate m-value between 'exact' and
                                        'mcmc' (default=exact)
    -mvalue_p_thres <FLOAT>             Compute m-values only for SNPs whose FE or RE2 p-value is
                                        below this threshold (default=1E-7)
    -mvalue_prior_alpha <ALPHA>         Alpha value for Beta dist prior Betadist(alpha,                                     beta) for existence of effect(default=1.0)
    -mvalue_prior_beta <BETA>           Beta value for Beta dist prior Betadist(alpha,                                      beta) for existence of effect (default=1.0,1.0)
    -mvalue_prior_sigma <FLOAT>         Sigma value for normal prior N(0, sigma^2) for effect
                                        (default=0.2)
    -mcmc_sample <INT>                  (MCMC) Number of samples (default=10,000)
    -mcmc_burnin <INT>                  (MCMC) Number of burn-in (default=1,000)
    -mcmc_max_num_flip <INT or FLOAT>   (MCMC) Usual move is flipping N bits where N ~
                                        U(1,max_num_flip). If an integer value i >= 1 is given,
                                        max_num_flip = i. If a float value 0 < k < 1 is given,
                                        max_num_flip = k * #studies. (default=0.1)
    -mcmc_prob_random_move <FLOAT>      (MCMC) Probability that a complete randomization move is
                                        suggested (default=0.01)
    -seed <INT>                         Random number generator seed (default=0)
    -thread <INT>                       Number of thread to use(default=1)
    -verbose                            Print RSID verbosely per every 1,000 SNPs (default=false)
    -help                               Print help
```
#
## Prerequisite<a name="prerequisite"></a>
In order to Improve and Edit METASOFT C++, <a href="https://www.boost.org/">boost library</a> is required.
```
boost-x64 1.71.0
```

#
## Reference <a name="reference"></a>

<a href="http://genetics.cs.ucla.edu/meta/">Original Software </a>