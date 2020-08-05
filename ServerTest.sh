echo "thread10"
date
sudo ./MEGASOFT/MEGASOFT.o --input ./TestData/performance_test/inputMS_1000x40.txt --pvalue_table ./TestData/HanEskinPvalueTable.txt --output ./ServerTestResult/posterior_1000x40_10.txt --log ./TestData/log_c++_static.txt --mvalue --mvalue_method mcmc --mcmc_sample 1000000 --seed 0 --mvalue_prior_sigma 0.05 --mvalue_prior_alpha 1 --mvalue_prior_beta 5 --mvalue_p_thres 1.0 --thread 10 > /dev/null
date
echo "thread20"
date
sudo ./MEGASOFT/MEGASOFT.o --input ./TestData/performance_test/inputMS_1000x40.txt --pvalue_table ./TestData/HanEskinPvalueTable.txt --output ./ServerTestResult/posterior_1000x40_20.txt --log ./TestData/log_c++_static.txt --mvalue --mvalue_method mcmc --mcmc_sample 1000000 --seed 0 --mvalue_prior_sigma 0.05 --mvalue_prior_alpha 1 --mvalue_prior_beta 5 --mvalue_p_thres 1.0 --thread 20 > /dev/null
date
echo "thread30"
date
sudo ./MEGASOFT/MEGASOFT.o --input ./TestData/performance_test/inputMS_1000x40.txt --pvalue_table ./TestData/HanEskinPvalueTable.txt --output ./ServerTestResult/posterior_1000x40_30.txt --log ./TestData/log_c++_static.txt --mvalue --mvalue_method mcmc --mcmc_sample 1000000 --seed 0 --mvalue_prior_sigma 0.05 --mvalue_prior_alpha 1 --mvalue_prior_beta 5 --mvalue_p_thres 1.0 --thread 30 >> /dev/null
date