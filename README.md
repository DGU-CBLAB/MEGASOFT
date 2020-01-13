# MEGASOFT (c++)

## Table of Content
- [About](#about)
- [How to Use](#use)
- [Prerequisite](#prerequisite)
- [Compile](#compile) ([windows](#windows), [linux](#linux))
- [Reference](#reference)

#
## About <a name="about"></a>
<a href="https://github.com/JuhunC/MEGASOFT">MEGASOFT</a> is a next-level software of <a href="http://genetics.cs.ucla.edu/meta/">METASOFT</a>. In Genome-wide association Study Analysis(GWAS), METASOFT is considered a widely used software, which is effective in performing wide range of basic and complicated meta-analytic methods. However, in terms of scalability, METASOFT lack of memory efficiency as well as computing performance. <a href="https://github.com/JuhunC/MEGASOFT/">MEGASOFT</a> is an improved version of METASOFT refered as <a href="https://github.com/JuhunC/MEGASOFT">MEGASOFT</a> where load balancing is used to improve computing performance as well as memory efficiency.

#

## How to Use MEGASOFT<a name="use"></a>

METASOFT <a href="http://genetics.cs.ucla.edu/meta/">Overview Link</a>

User's MEGASOFT Guide
```
MEGASOFT_WIN64.exe -input [inputMS.txt] [option1] [option2] ...
MEGASOFT_LINUX.o -input [inputMS.txt] [option1] [option2] ...
```
```
-input <FILE>                       Input file (Required)
```
```
-output <FILE>                      Output file (default='out')
```
```
-log <FILE>                         Log file (default='log')
```
```
-pvalue_table <FILE>                Pvalue table file (default='HanEskinPvalueTable.txt')
```
```
-lambda_mean <FLOAT>                (Random Effects) User-specified lambda for mean effect part(default=1.0)
```
```
-lambda_hetero <FLOAT>              (Random Effects) User-specified lambda for heterogeneity part (default=1.0)
```
```
-binary_effects                     Compute binary effects model p-value (default=false)
```
```
-binary_effects_sample <INT>        (Binary effects) Number of importance sampling samples(default=1,000)
```
```
-binary_effects_p_thres <FLOAT>     (Binary effects) P-value threshold determining if we will use large number of samples (default=1E-4)
```
```
-binary_effects_large <INT>         (Binary effects) Large number of importance sampling samples for p-values above threshold (default=100,000)
```
```
-mvalue                             Compute m-value (default=false)
```
```
-mvalue_method <METHOD_NAME>        Which method to use to calculate m-value between 'exact' and 'mcmc' (default=exact)
```
```
-mvalue_p_thres <FLOAT>             Compute m-values only for SNPs whose FE or RE2 p-value is below this threshold (default=1E-7)
```
```
-mvalue_prior_alpha <ALPHA>         Alpha value for Beta dist  prior Betadist(alpha,beta) for existence of effect (default=1.0)
```
```
-mvalue_prior_beta <BETA>           Beta value for Beta dist  prior Betadist(alpha,beta) for existence of effect(default=1.0,1.0)
```
```
-mvalue_prior_sigma <FLOAT>         Sigma value for normal  prior N(0, sigma^2) for effect(default=0.2)
```
```
-mcmc_sample <INT>                  (MCMC) Number of samples  (default=10,000)
```
```
-mcmc_burnin <INT>                  (MCMC) Number of burn-in  (default=1,000)
```
```
-mcmc_max_num_flip <INT or FLOAT>   (MCMC) Usual move is  flipping N bits where N ~ U(1,max_num_flip). If an integer value i >= 1 is given, max_num_flip = i. If a float value 0 < k < 1 is given, max_num_flip = k * #studies. (default=0.1)
```
```
-mcmc_prob_random_move <FLOAT>      (MCMC) Probability that a complete randomization move is suggested (default=0.01)
```
```
-seed <INT>                         Random number generator seed (default=0)
```
```
-thread <INT>                       Number of thread to use(default=1)
```
```
-verbose                            Print RSID verbosely per every 1,000 SNPs (default=false)
```
```
-help                               Print help
```
#
## Prerequisite<a name="prerequisite"></a>
In order to Improve and Edit MEGASOFT, <a href="https://www.boost.org/">boost library</a> is required.
```
Windows                             boost-x64 1.72.0
Linux                               boost-x64 1.72.0
```
#
## How to Compile <a href="https://github.com/JuhunC/MEGASOFT">MEGASOFT</a><a name="compile"></a>
<a href="https://www.boost.org/users/download/">Boost Download link</a>

### Windows x64<a name="windows"></a> - using Visual Studio 2019(v142)
1. `boostrap.bat`

2. `b2 variant=debug,release link=static threading=multi address-model=64 runtime-link=static -j4 install --prefix=stage`
3. `Project Properties Page -> Configuration Properties -> VC++ Directories -> Include Directories`
    ``` 
    ADD YOUR_BOOST_DIR 
    ..\boost_1_72_0\stage\include\boost-1_72\
    ```
4. `Project Properties Page -> Configuration Properties -> C/C++ -> Code Generation -> Runtime Library`
    ```
    Set to "Multi-threaded(/MT)"
    ```
5. `Project Properties Page -> Configuration Properties -> Linker -> General -> Additional Library Directories`
    ```
    ADD 
    YOUR_BOOST_LIB_DIR ..\boost_1_72_0\stage\lib\
    ```
6. `Project Properties Page -> Configuration Properties -> Linker -> Input -> Additional Dependencies`
    ```
    ADD
    libboost_program_options-vc142-mt-s-x64-1_72.lib
    libboost_thread-vc142-mt-s-x64-1_72.lib
    libboost_system-vc142-mt-s-x64-1_72.lib
    ```
### Linux - using MinGW(g++)<a name="linux"></a>
1. `boostrap.sh`

2. `b2 link=static`

3.  ```
    export BOOST_ROOT=/home/user/Desktop/boost_1_72_0
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BOOST_ROOT/stage/lib
    ```
4. Final Compile using g++
    ```
    sudo g++ -pthread -I BOOST_ROOT MetaSnp.cpp  MetaSnp.h MetaSoft.cpp LD_LIBRARY_PATH/libboost_thread.a LD_LIBRARY_PATH/libboost_program_options.a LD_LIBRARY_PATH/libboost_system.a -o out.o -g3 -static
    ```


## Reference <a name="reference"></a>

- <a href="http://genetics.cs.ucla.edu/meta/">METASOFT</a>
- <a href="https://www.boost.org/">Boost</a>
