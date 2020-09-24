#include "MetaSnp.h"
#include <ctime>
#include <boost/program_options.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include<mutex>
#include<condition_variable>
namespace po = boost::program_options;

struct{
	std::vector<std::string> lines;
	FILE* outFile;
	unsigned int* done_count_ptr;
}typedef thread_struct;

// Multi-Thread variables
static int threadNum_ = 1;

static boost::mutex mtx;
static boost::mutex done_mtx;
static bool flag = false;

// Arguments and default values
static std::string  inputFile_ = "";
static std::string  outputFile_ = "out";
static std::string  logFile_ = ".\\log";
static std::string  pvalueTableFile_ = "HanEskinPvalueTable.txt";
static double  inputLambdaMeanEffect_ = 1.0;
static double  inputLambdaHeterogeneity_ = 1.0;
static bool willComputeMvalue_ = false;
static double  priorSigma_ = 0.2;
static double  priorAlpha_ = 1.0;
static double  priorBeta_ = 1.0;
static double  mvaluePvalueThreshold_ = 1E-7;
static std::string  _mvalueMethodmvalueMethod_ = "exact";
static long    mcmcSample_ = 10000;
static long    mcmcBurnin_ = 1000;
static double  mcmcProbRandom_ = 0.01;
static double  mcmcMaxNumFlip_ = 0.1;
static bool willComputeBinaryEffects_ = false;
static long    binaryEffectsSample_ = 1000;
static long    binaryEffectsLargeSample_ = 100000;
static double  binaryEffectsPvalueThreshold_ = 1E-4;
static unsigned int     seed_ = 0;
static bool isVerbose_ = false;

// Internally used variables
static int    numSnps_;
static int    maxNumStudy_;
static double outputLambdaMeanEffect_;
static double outputLambdaHeterogeneity_;
static std::vector<double> meanEffectParts_;
static std::vector<double> heterogeneityParts_;
static std::string argsSummary_;

const double expectedMedianHanEskinHeterogeneityPart_[] = // from nStudy 2 to 50
{ 0.2195907137,0.2471516439,0.2642270318,0.2780769264,0.2886280267,0.2977812664,0.3020051148,0.3091428179,0.3158605559,0.3221788173,0.3259133140,0.3295976587,0.3335375196,0.3358395088,0.3368309971,0.3421941686,0.3448030927,0.3463590948,0.3477384754,0.3487900288,0.3494938171,0.3542087791,0.3573286353,0.3589703411,0.3586951356,0.3596101209,0.3605611682,0.3624799993,0.3648322669,0.3659817739,0.3671267389,0.3693952373,0.3693395144,0.3696863113,0.3706067524,0.3718103285,0.3749536619,0.3758886239,0.3753612342,0.3781458299,0.3798346038,0.3763434983,0.3796968747,0.3784334922,0.3794411347,0.3808582942,0.3813485882,0.3843230993,0.3824863479 };

void printErrorAndQuit(std::string msg);

void handleArguments(int argc, char* argv[]) 
{
	boost::program_options::options_description desc("Allowed Options");
	desc.add_options()
		("help", "Print help")
		("input", po::value<std::string>(), "Input file (Required)")
		("output", po::value<std::string>(), "Output file (default='out')")
		("pvalue_table", po::value<std::string>(), "value table file (default='HanEskinPvalueTable.txt')")
		("log", po::value<std::string>(), "Log file (default='log')")
		("lambda_mean", po::value<double>(), "(Random Effects) User-specified lambda for mean effect part (default=1.0)")
		("lambda_hetero", po::value<double>(), "(Random Effects) User-specified lambda for heterogeneity part (default=1.0)")
		("mvalue", po::bool_switch()->default_value(false), "Compute m-value(default=false)")
		("mvalue_prior_sigma", po::value<double>(), "Sigma value for normal prior N(0, sigma^2) for effect (default=0.2)")
		("mvalue_prior_alpha", po::value<std::string>(), "Alpha value for Beta dist prior Betadist(alpha,beta) for existence of effect (default=1.0,1.0)")
		("mvalue_prior_beta", po::value<std::string>(), "Beta value for Beta dist prior Betadist(alpha,beta) for existence of effect (default=1.0,1.0)")
		("mvalue_p_thres", po::value<double>(), "Compute m-values only for SNPs whose FE or RE2 p-value is below this threshold (default=1E-7)")
		("mvalue_method", po::value<std::string>(), "Which method to use to calculate m-value between 'exact' and 'mcmc' (default=exact)")
		("mcmc_sample", po::value<long>(), "(MCMC) Number of samples (default=10,000)")
		("mcmc_burnin", po::value<long>(), "(MCMC) Number of burn-in (default=1,000)")
		("mcmc_prob_random_move", po::value<double>(), "(MCMC) Probability that a complete randomization move is suggested (default=0.01)")
		("mcmc_max_num_flip", po::value<double>(), "(MCMC) Usual move is flipping N bits where N ~ U(1,max_num_flip). If an integer value i >= 1 is given, max_num_flip = i. If a float value 0 < k < 1 is given, max_num_flip = k * #studies. (default=0.1)")
		("binary_effects", po::value<std::string>(), "Compute binary effects model p-value (default=false)")
		("binary_effects_sample", po::value<long>(), "(Binary effects) Number of importance sampling samples (default=1,000)")
		("binary_effects_large", po::value<long>(), "(Binary effects) Large number of importance sampling samples for p-values above threshold (default=100,000)")
		("binary_effect_p_thres", po::value<double>(), "(Binary effects) P-value threshold determining if we will use large number of samples (default=1E-4)")
		("seed", po::value<int>(), "Random number generator seed (default=0)")
		("verbose", po::bool_switch()->default_value(false), "Print RSID verbosely per every 1,000 SNPs (default=false)")
		("thread", po::value<int>(), "Set Number of Threads")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	

	if (vm.count("help")) 
	{
		std::cout << desc << std::endl;
		std::exit(1);
	}
	if (vm.count("input")) 
	{
		inputFile_ = vm["input"].as<std::string>();
	}
	
	if (vm.count("output")) 
	{
		outputFile_ = vm["output"].as<std::string>();
	}
	if (vm.count("pvalue_table")) 
	{
		pvalueTableFile_ = vm["pvalue_table"].as<std::string>();
	}
	if (vm.count("log")) 
	{
		logFile_ = vm["log"].as<std::string>();
	}
	if (vm.count("lambda_mean")) 
	{
		inputLambdaMeanEffect_ = vm["lambda_mean"].as<double>();
	}
	if (vm.count("lambda_hetero")) 
	{
		inputLambdaHeterogeneity_ = vm["lambda_hetero"].as<double>();
	}
	if (vm.count("mvalue")) 
	{
		willComputeMvalue_ = true;
		if (vm.count("mvalue_prior_sigma")) 
		{
			priorSigma_ = vm["mvalue_prior_sigma"].as<double>();
		}
		if (vm.count("mvalue_prior_alpha")) 
		{
			priorAlpha_ = stod(vm["mvalue_prior_alpha"].as<std::string>());
		}
		if (vm.count("mvalue_prior_beta")) 
		{
			priorBeta_ = stod(vm["mvalue_prior_beta"].as<std::string>());
		}
		if (vm.count("mvalue_p_thres")) 
		{
			mvaluePvalueThreshold_ = vm["mvalue_p_thres"].as<double>();
		}
		if (vm.count("mvalue_method")) 
		{
			mvalueMethod_ = vm["mvalue_method"].as<std::string>();
		}
		if (mvalueMethod_ == "mcmc") 
		{
			if (vm.count("mcmc_sample")) 
			{
				mcmcSample_ = vm["mcmc_sample"].as<long>();
			}
			if (vm.count("mcmc_burnin")) 
			{
				mcmcBurnin_ = vm["mcmc_burnin"].as<long>();
			}
			if (vm.count("mcmc_prob_random_move")) 
			{
				mcmcProbRandom_ = vm["mcmc_prob_random_move"].as<double>();
			}
			if (vm.count("mcmc_max_num_flip")) 
			{
				mcmcMaxNumFlip_ = vm["mcmc_max_num_flip"].as<double>();
			}
		}
	}
	if (vm.count("binary_effects")) 
	{
		willComputeBinaryEffects_ = true;
		if (vm.count("binary_effects_sample")) 
		{
			binaryEffectsSample_ = vm["binary_effects_sample"].as<long>();
		}
		if (vm.count("binary_effects_large")) 
		{
			binaryEffectsLargeSample_ = vm["binary_effects_large"].as<long>();
		}
		if (vm.count("binary_effects_p_thres")) 
		{
			binaryEffectsPvalueThreshold_ = vm["binary_effects_p_thres"].as<double>();
		}
	}
	if (vm.count("seed")) 
	{
		seed_ = vm["seed"].as<int>();
	}
	if (vm["verbose"].as<bool>() == true) 
	{
		isVerbose_ = true;
	}
	if (vm.count("thread")) 
	{
		threadNum_ = vm["thread"].as<int>();
#ifdef FORCE_THREAD
		threadNum_ = THREAD;
#endif
	}
	if (vm.count("help")) 
	{
		std::cout << "------------------------------------------------"<< std::endl;
		std::cout << "The format of input_file:" << std::endl;
		std::cout << "  Each row is each SNP." << std::endl;
		std::cout << "  1st column is RSID." << std::endl;
		std::cout << "  2nd and 3rd column are beta and its standard error for 1st study." << std::endl;
		std::cout << "  4th and 5th column are beta and its standard error for 2nd study." << std::endl;
		std::cout << "  6th and 7th column are beta and its standard error for 3rd study." << std::endl;
		std::cout << "  and so on..." << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		std::cout << std::endl;
		std::exit(-1);
	}


	std::cout << "-------- Processing arguments ---------" << std::endl;
	if (inputFile_ == "") 
	{
		printErrorAndQuit("A valid input file must be specified using option -input");
	}
	if (inputLambdaMeanEffect_ <= 0.0) 
	{
		printErrorAndQuit("lambda_mean option takes a float value > 0");
	}
	if (inputLambdaHeterogeneity_ <= 0.0) 
	{
		printErrorAndQuit("lambda_hetero option takes a float value > 0");
	}
	if (priorSigma_ <= 0.0) 
	{
		printErrorAndQuit("mvalue_prior_sigma option takes a float value > 0");
	}
	if (priorAlpha_ <= 0.0 || priorBeta_ <= 0.0) 
	{
		printErrorAndQuit("mvalue_prior_beta option takes two float values > 0");
	}
	if (mvaluePvalueThreshold_ < 0.0 || mvaluePvalueThreshold_ > 1.0) 
	{
		printErrorAndQuit("mvalue_p_thres takes a float value between 0 and 1");
	}
	if (mvalueMethod_ != "exact" &&
		mvalueMethod_ != "mcmc" &&
		mvalueMethod_ != "variational") 
	{
		printErrorAndQuit("mvalue_method option only takes a value 'exact' or 'mcmc'");
	}
	if (mcmcSample_ < 1) 
	{
		printErrorAndQuit("mcmc_sample option takes an integer value > 0");
	}
	if (mcmcBurnin_ < 1) 
	{
		printErrorAndQuit("mcmc_burnin option takes an integer value > 0");
	}
	if (mcmcSample_ < mcmcBurnin_) 
	{
		printErrorAndQuit("mcmc_sample must be larger than mcmc_burnin");
	}
	if (mcmcProbRandom_ < 0.0 || mcmcProbRandom_ > 1.0) 
	{
		printErrorAndQuit("mcmc_prob_random takes a float value between 0 and 1");
	}
	if (mcmcMaxNumFlip_ <= 0.0) 
	{
		printErrorAndQuit("mcmc_max_num_flip takes a value > 0");
	}
	if (binaryEffectsSample_ < 1) 
	{
		printErrorAndQuit("binary_effects_sample option takes an integer value > 0");
	}
	if (binaryEffectsLargeSample_ < 1) 
	{
		printErrorAndQuit("binary_effects_large option takes an integer value > 0");
	}
	if (binaryEffectsPvalueThreshold_ < 0.0 || binaryEffectsPvalueThreshold_ > 1.0) 
	{
		printErrorAndQuit("binary_effects_p_thres takes a float value between 0 and 1");
	}
	
	// Make summary for printing
	
	argsSummary_ = "";
	for (int i = 0; i < argc;i++) 
	{
		argsSummary_ += std::string(argv[i]) + " ";
	}
	argsSummary_ += "\n";
}

void printErrorAndQuit(std::string msg) 
{
	printf("%s\n",msg.c_str());
	std::exit(-1);
}
void* thr_func(void* args)
{
	thread_struct* thr_args = (thread_struct*)args;
	std::vector<std::string> LINES = thr_args->lines;
	FILE* outFile = thr_args->outFile;
	MetaSnp* metaSnp;    // Store only 1 Snp at a time in memory.

	for (std::string line : LINES)
	{
		std::vector<std::string> tokens;
		split(tokens, line, " ");
		if (tokens.size() > 1) 
		{             // only if non-empty
			if (tokens.at(0).at(0) != '#') 
			{ // only if non-comment
				std::string rsid = tokens.at(0);
				metaSnp = new MetaSnp(rsid);
				if (tokens.size() % 2 == 0)
				{
					printf("WARNING: # of Columns must be odd including Rsid. Last column is ignored.");
				}

				int nStudy = (int)((tokens.size() - 1) / 2);
				if (nStudy > maxNumStudy_) 
				{
					maxNumStudy_ = nStudy;
				}

				for (int i = 0; i < nStudy; i++) 
				{
					double beta;
					double standardError;
					if (tokens.at((2.0 * i) + 1.0).compare("NA") == 0 ||
						tokens.at((2.0 * i) + 1.0).compare("N/A") == 0 ||
						tokens.at((2.0 * i) + 1.0).compare("NA") == 0 ||
						tokens.at((2.0 * i) + 1.0).compare("N/A") == 0) 
					{
						metaSnp->addNaStudy();
					}
					else 
					{
						try {
							beta = stod(tokens.at((2.0 * i) + 1.0));
							standardError = stod(tokens.at(2.0 * i + 2.0));
							if (standardError <= 0.0) 
							{
								printf("Standard error cannot be <= zero (%f th column is %f) in the following line.\n",
									2.0 * i + 3.0, standardError);
								printf("%s", line.c_str());
								std::exit(-1);
							}
							metaSnp->addStudy(beta, standardError);
						}
						catch (std::exception es) 
						{
							printf("Incorrect float value in following line. Possibly not a double");
							printf("%s", line.c_str());
							std::exit(-1);
						}
					}
				}
				if (metaSnp->getNStudy() > 1) 
				{
					metaSnp->initMem();
					// Analyze 1 Snp on-the-fly.
					if (isVerbose_ && numSnps_ % 1000 == 0) 
					{
						printf("Analyzing SNP #%d (%s)\n", numSnps_ + 1, rsid.c_str());
					}
					// FE, RE, and New RE
					metaSnp->computeFixedEffects();
					metaSnp->computeRandomEffects();
					metaSnp->computeHanEskin(inputLambdaMeanEffect_, inputLambdaHeterogeneity_);

					meanEffectParts_.push_back(metaSnp->getStatisticHanEskinMeanEffectPart());
					
					double h = metaSnp->getStatisticHanEskinHeterogeneityPart();
					
					if (h > 0.0) 
					{
						heterogeneityParts_.push_back(h);
					}
					// Binary effects model
					if (willComputeBinaryEffects_) 
					{
						metaSnp->computeBinaryEffectsPvalue(binaryEffectsSample_, rand());
						if (metaSnp->getPvalueBinaryEffects() <= binaryEffectsPvalueThreshold_) 
						{
							metaSnp->computeBinaryEffectsPvalue(binaryEffectsLargeSample_, rand());
						}
					}
					// Mvalues
					if (willComputeMvalue_) 
					{
						if (metaSnp->getPvalueFixedEffects() <= mvaluePvalueThreshold_ ||
							metaSnp->getPvalueHanEskin() <= mvaluePvalueThreshold_) 
						{
							if (mvalueMethod_.compare("exact") == 0) 
							{
								metaSnp->computeMvalues(priorAlpha_, priorBeta_, priorSigma_);
							}
							else if (mvalueMethod_.compare("mcmc") == 0) 
							{
								metaSnp->computeMvaluesMCMC(priorAlpha_, priorBeta_, priorSigma_,
									mcmcSample_, mcmcBurnin_, mcmcProbRandom_,
									mcmcMaxNumFlip_,
									rand());
							}
							else 
							{
								std::cout << mvalueMethod_ << std::endl;
								assert(false);
							}
						}
					}
					numSnps_++;
				}// end of if
				mtx.lock();
				metaSnp->printResults(outFile);
				mtx.unlock();
			}// end of if('#')
		}// end of if(token.size()>1)
		tokens.clear();

		done_mtx.lock();
		*thr_args->done_count_ptr = *thr_args->done_count_ptr+1;
		done_mtx.unlock();
	}
}
void doMetaAnalysis() 
{
	srand(seed_);
	numSnps_ = 0;
	maxNumStudy_ = 0;
	meanEffectParts_ = std::vector<double>();
	heterogeneityParts_ = std::vector<double>();
	FILE* inFile, *outFile;
	try {
		inFile = fopen(inputFile_.c_str(), "r");
	}
	catch (std::exception e) 
	{
		printErrorAndQuit("ERROR: Input file cannot be opened");
	}
	try {
		outFile = fopen(outputFile_.c_str(), "w");
	}
	catch (std::exception e) 
	{
		printErrorAndQuit("ERROR: Ouput file cannot be opened");
	}

	// Print headings
	MetaSnp::printHeadings(outFile);

	try {
		std::ifstream inStream(inputFile_);
		std::vector<std::string> LINES;
		unsigned int done_count = 0;
		std::string tmp_ln;
		thread_struct* args;

#ifdef LINUX
		int err;
		pthread_t* thr_ptr = new pthread_t();
		std::vector<pthread_t*> thread_pool;
#elif WINDOWS
		std::thread* thr_ptr;
		std::vector<std::thread*> thread_pool;
#endif

		while (std::getline(inStream, tmp_ln))
		{
			LINES.push_back(tmp_ln);
		}

		int N = LINES.size() / threadNum_ + (LINES.size()%threadNum_ != 0);

		//TODO: Divide Lines;
		for (int nidx = 0; nidx < threadNum_; nidx++)
		{
			std::vector<std::string> args_lines;
			for (int i = N * nidx; i < N * (nidx + 1) && i < LINES.size(); i++) 
			{
				args_lines.push_back(LINES.at(i));
			}
			args = new thread_struct();
			args->outFile = outFile;
			args->done_count_ptr = &done_count;
			args->lines = args_lines;;

		#ifdef LINUX
			err = pthread_create(thr_ptr, NULL, &thr_func, args);
			if (err != 0)
			{
				// printf("\nCan't create thread: [%s]", strerr(err));
				exit(ERROR_THREAD_CREATE);
			}
		#elif WINDOWS
			thr_ptr = new std::thread(thr_func, args);
		#endif

			thread_pool.push_back(thr_ptr);
		}

		while (done_count < LINES.size()) {
			//std::this_thread::sleep_for(std::chrono::microseconds(100));
			std::cout << "Current Progress : " << done_count <<"/"<<LINES.size()<< " finished." << "\r";
		}




		// Join Threads
		for (int nidx = 0; nidx < thread_pool.size(); nidx++)
		{
		#ifdef LINUX
			err = pthread_join(*(pthread_t*)thread_pool.at(nidx), NULL);
			if (err != 0)
			{
				// printf("\nCan't join thread: [%s]\n", strerror(err));
				exit(ERROR_THREAD_JOIN);
			}
		#elif WINDOWS
			thread_pool.at(nidx)->join();
		#endif
		}// end of for(nidx) : Thread Join


		try {
			fclose(inFile);
			fclose(outFile);
		}
		catch (std::exception e)
		{
			printf("ERROR: file cannot be closed");
			std::exit(ERROR_IO_FILE_CLOSE);
		}
	}
	catch (std::exception e) 
	{
		printf("ERROR: doMetaAnalysis [%s]", e.what());
		std::exit(ERROR_META_ANALYSIS);
	}
	

	// Reorder the Result File
	try {
		std::string tmp_str;
		FILE* file = fopen(outputFile_.c_str(), "r");
		std::ifstream Instream(outputFile_.c_str());

		std::getline(Instream, tmp_str); // ignore first line
		
		std::vector<ThreadResult> total;

		while (std::getline(Instream, tmp_str)) 
		{
			ThreadResult mt(stoi(tmp_str.substr(0, tmp_str.find('\t'))), tmp_str);
			total.push_back(mt);
		}

		sort(total.begin(), total.end(), compareKeyValues);
		Instream.close();
		FILE* outfile = fopen(outputFile_.c_str(), "w");
		MetaSnp::printHeadings(outfile);
		for (int i = 0; i < total.size(); i++) 
		{
			fprintf(outfile, "%s\n", total.at(i).values_str.c_str());
		}
		
		fclose(outfile);
	}
	catch (std::exception e) 
	{
		printf("ERROR: Posterior.txt file has been altered!\n[%s]",e.what());
		std::exit(-1);
	}
}

void computeLambda() 
{
	double median;
	double expectedMedian;
	if (meanEffectParts_.size() > 0) 
	{
		std::sort(meanEffectParts_.begin(), meanEffectParts_.end());
		median = meanEffectParts_.at((int)(meanEffectParts_.size() / 2.0));
		expectedMedian = pow(boost::math::find_location<boost::math::normal>(boost::math::complement(0, 0.25, 1.0)), 2.0);
		outputLambdaMeanEffect_ = median / expectedMedian;
	}
	if (heterogeneityParts_.size() > 0) 
	{
		std::sort(heterogeneityParts_.begin(), heterogeneityParts_.end());
		median = heterogeneityParts_.at((int)(heterogeneityParts_.size() / 2.0));
		if (maxNumStudy_ > 50) 
		{
			expectedMedian = pow(boost::math::find_location<boost::math::normal>(boost::math::complement(0, 0.25, 1.0)), 2.0);
		}
		else 
		{
			expectedMedian = expectedMedianHanEskinHeterogeneityPart_[maxNumStudy_ - 2];
		}
		outputLambdaHeterogeneity_ = median / expectedMedian;
	}
}

void printLog(double time) 
{
	try {
		FILE* outFile = fopen(logFile_.c_str(),"w");
		
		fprintf(outFile,"Arguments: %s ", argsSummary_.c_str());
		fprintf(outFile,"Input File: %s\n", inputFile_.c_str());
		fprintf(outFile,"Output File: %s\n", outputFile_.c_str());
		fprintf(outFile,"Log File: %s\n", logFile_.c_str());
		fprintf(outFile,"p-value Table File: %s\n", pvalueTableFile_.c_str());
		fprintf(outFile,"Number of SNPs analyzed: %d\n", numSnps_);
		fprintf(outFile,"Maximum number of studies: %d\n", maxNumStudy_);
		fprintf(outFile,"Specified lambda for   mean effect part (default = 1.0): %f\n", inputLambdaMeanEffect_);
		fprintf(outFile,"Specified lambda for heterogeneity part (default = 1.0): %f\n", inputLambdaHeterogeneity_);
		fprintf(outFile,"Newly calculated inflation factor lambda for   mean effect part: %f\n", outputLambdaMeanEffect_);
		fprintf(outFile,"Newly calculated inflation factor lambda for heterogeneity part: %f\n", outputLambdaHeterogeneity_);
		fprintf(outFile,"Elapsed time : %.3f minutes\n",time);
		fclose(outFile);
	}
	catch (std::exception e) 
	{
		printf("ERROR: error encountered while writing in log file");
		std::exit(-1);
	}
}


int main(int argc, char* argv[]) 
{	
	time_t startTime = time(NULL);

	handleArguments(argc, argv);
	std::cout << "Arguments: " + argsSummary_ << std::endl;

	MetaSnp::readPvalueTableFile(pvalueTableFile_);

	std::cout<<"----- Performing meta-analysis\n";
	doMetaAnalysis();

	std::cout << "---- Performing lambda compute\n";
	computeLambda();

	std::cout << "---- print Log\n";
	time_t endTime = time(NULL);
	printLog(difftime(endTime, startTime)/60.0);

	std::cout<<"----- Finished\n";
	
	std::cout << "----- Elapsed time: " << difftime(endTime,startTime)/60.0 << " minutes\n";

	return DONE_NORMAL;
}
