echo "thread1"
date
sudo ./MEGASOFT/MEGASOFT.o --input ./TestData/performance_test/inputMS_1000x40.txt --pvalue_table ./TestData/HanEskinPvalueTable.txt --output ./ServerTestResult/posterior_1000x40_1.txt --log ./TestData/log_c++_static.txt --mvalue --mvalue_method mcmc --mcmc_sample 1000000 --seed 0 --mvalue_prior_sigma 0.05 --mvalue_prior_alpha 1 --mvalue_prior_beta 5 --mvalue_p_thres 1.0 --thread 1 > /dev/null
date
echo "thread3"
date
sudo ./MEGASOFT/MEGASOFT.o --input ./TestData/performance_test/inputMS_1000x40.txt --pvalue_table ./TestData/HanEskinPvalueTable.txt --output ./ServerTestResult/posterior_1000x40_3.txt --log ./TestData/log_c++_static.txt --mvalue --mvalue_method mcmc --mcmc_sample 1000000 --seed 0 --mvalue_prior_sigma 0.05 --mvalue_prior_alpha 1 --mvalue_prior_beta 5 --mvalue_p_thres 1.0 --thread 3 > /dev/null
date
echo "thread5"
date
sudo ./MEGASOFT/MEGASOFT.o --input ./TestData/performance_test/inputMS_1000x40.txt --pvalue_table ./TestData/HanEskinPvalueTable.txt --output ./ServerTestResult/posterior_1000x40_5.txt --log ./TestData/log_c++_static.txt --mvalue --mvalue_method mcmc --mcmc_sample 1000000 --seed 0 --mvalue_prior_sigma 0.05 --mvalue_prior_alpha 1 --mvalue_prior_beta 5 --mvalue_p_thres 1.0 --thread 5 >> /dev/null
date