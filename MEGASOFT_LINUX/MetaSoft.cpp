#include "MetaSnp.h"
#include <ctime>
#include <boost/program_options.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include<mutex>
#include<condition_variable>

namespace po = boost::program_options;

typedef struct
{
	FILE* outFile;
	std::vector<std::string> lines;
	unsigned int* done_count_ptr;
}thread_input;

static unsigned int _threadNum = 1;

static boost::mutex 	Mtx;
static boost::mutex		DoneMtx;

static std::string	_inputFile 						= "";
static std::string 	_outputFile						= "out";
static std::string 	_logFile 						= ".\\log";
static std::string	_pvalueTableFile 				= "HanEskinPvalueTable.txt";

static float 		_inputLambdaMeanEffects			= 1.0;
static float		_inputLambdaHeterogeneity		= 1.0;
static float 		_priorSigma						= 0.2;
static float		_priorAlpha						= 1.0;
static float 		_priorBeta						= 1.0;
static float 		_mvaluePvalueThreshold			= 1E-7;
static std::string	_mvalueMethod					= "exact";
static unsigned int	_mcmcSample						= 10_000;
static unsigned int	_mcmcBurnin						= 1_000;
static float 		_mcmcProbRandom					= 0.01;
static float 		_mcmcMaxNumFlip					= 0.1;
static unsigned int	_binaryEffectsSample			= 1_000;
static unsigned int	_binaryEffectsLargeSample		= 100_000;
static float 		_binaryEffectsPvalueThreshold	= 1E-4;

static bool			_willComputeMvalue 				= false;
static bool			_willComputeBinaryEffects		= false;

static unsigned int _seed							= 0;
static bool 		_verbose						= false;

static unsigned int	_numSNPs;
static unsigned int	_maxNumStudy;
static float 		_outputLambdaMeanEffects;
static float 		_outputLambdaHeterogeneity;
static std::vector<float>	_meanEffectsParts;
static std::vector<float> 	_heterogeneityParts;
static std::string			_argsSummary;

const float expectedMedianHanEskinHeterogeneityPart_[] = // from nStudy 2 to 50
{ 0.2195907137,0.2471516439,0.2642270318,0.2780769264,0.2886280267,0.2977812664,0.3020051148,0.3091428179,0.3158605559,0.3221788173,0.3259133140,0.3295976587,0.3335375196,0.3358395088,0.3368309971,0.3421941686,0.3448030927,0.3463590948,0.3477384754,0.3487900288,0.3494938171,0.3542087791,0.3573286353,0.3589703411,0.3586951356,0.3596101209,0.3605611682,0.3624799993,0.3648322669,0.3659817739,0.3671267389,0.3693952373,0.3693395144,0.3696863113,0.3706067524,0.3718103285,0.3749536619,0.3758886239,0.3753612342,0.3781458299,0.3798346038,0.3763434983,0.3796968747,0.3784334922,0.3794411347,0.3808582942,0.3813485882,0.3843230993,0.3824863479 };

void handleArguments(int argc, char* argv[]);
void errorHandler(std::string msg);
void doMetaAnalysis();

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

//////////////////////////////////////////////////////////////////////////////////////////////////


void handleArguments(int argc, char* argv[]) 
{
	boost::program_options::options_description desc("Allowed Options");
	desc.add_options()
		("help", 					"Print help")
		("input", 					po::value<std::string>(), 					"Input file (Required)")
		("output", 					po::value<std::string>(), 					"Output file (default='out')")
		("pvalue_table", 			po::value<std::string>(), 					"value table file (default='HanEskinPvalueTable.txt')")
		("log", 					po::value<std::string>(), 					"Log file (default='log')")
		("lambda_mean", 			po::value<float>(), 						"(Random Effects) User-specified lambda for mean effect part (default=1.0)")
		("lambda_hetero", 			po::value<float>(), 						"(Random Effects) User-specified lambda for heterogeneity part (default=1.0)")
		("mvalue", 					po::bool_switch()->default_value(false), 	"Compute m-value(default=false)")
		("mvalue_prior_sigma", 		po::value<float>(), 						"Sigma value for normal prior N(0, sigma^2) for effect (default=0.2)")
		("mvalue_prior_alpha", 		po::value<std::string>(), 					"Alpha value for Beta dist prior Betadist(alpha,beta) for existence of effect (default=1.0,1.0)")
		("mvalue_prior_beta", 		po::value<std::string>(), 					"Beta value for Beta dist prior Betadist(alpha,beta) for existence of effect (default=1.0,1.0)")
		("mvalue_p_thres", 			po::value<float>(), 						"Compute m-values only for SNPs whose FE or RE2 p-value is below this threshold (default=1E-7)")
		("mvalue_method", 			po::value<std::string>(), 					"Which method to use to calculate m-value between 'exact' and 'mcmc' (default=exact)")
		("mcmc_sample", 			po::value<unsigned int>(), 					"(MCMC) Number of samples (default=10,000)")
		("mcmc_burnin", 			po::value<unsigned int>(), 					"(MCMC) Number of burn-in (default=1,000)")
		("mcmc_prob_random_move", 	po::value<float>(), 						"(MCMC) Probability that a complete randomization move is suggested (default=0.01)")
		("mcmc_max_num_flip", 		po::value<float>(), 						"(MCMC) Usual move is flipping N bits where N ~ U(1,max_num_flip). If an integer value i >= 1 is given, max_num_flip = i. If a float value 0 < k < 1 is given, max_num_flip = k * #studies. (default=0.1)")
		("binary_effects", 			po::value<std::string>(), 					"Compute binary effects model p-value (default=false)")
		("binary_effects_sample", 	po::value<unsigned int>(), 							"(Binary effects) Number of importance sampling samples (default=1,000)")
		("binary_effects_large", 	po::value<unsigned int>(), 							"(Binary effects) Large number of importance sampling samples for p-values above threshold (default=100,000)")
		("binary_effect_p_thres", 	po::value<float>(), 						"(Binary effects) P-value threshold determining if we will use large number of samples (default=1E-4)")
		("seed", 					po::value<unsigned int>(), 							"Random number generator seed (default=0)")
		("verbose", 				po::bool_switch()->default_value(false), 	"Print RSID verbosely per every 1,000 SNPs (default=false)")
		("thread", 					po::value<unsigned int>(), 							"Set Number of Threads")
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
		_inputFile = vm["input"].as<std::string>();
	}
	
	if (vm.count("output")) 
	{
		_outputFile = vm["output"].as<std::string>();
	}
	if (vm.count("pvalue_table")) 
	{
		_pvalueTableFile = vm["pvalue_table"].as<std::string>();
	}
	if (vm.count("log")) 
	{
		_logFile = vm["log"].as<std::string>();
	}
	if (vm.count("lambda_mean")) 
	{
		_inputLambdaMeanEffects = vm["lambda_mean"].as<float>();
	}
	if (vm.count("lambda_hetero")) 
	{
		_inputLambdaHeterogeneity = vm["lambda_hetero"].as<float>();
	}
	if (vm.count("mvalue")) 
	{
		_willComputeMvalue = true;
		if (vm.count("mvalue_prior_sigma")) 
		{
			_priorSigma = vm["mvalue_prior_sigma"].as<float>();
		}
		if (vm.count("mvalue_prior_alpha")) 
		{
			_priorAlpha = stod(vm["mvalue_prior_alpha"].as<std::string>());
		}
		if (vm.count("mvalue_prior_beta")) 
		{
			_priorBeta = stod(vm["mvalue_prior_beta"].as<std::string>());
		}
		if (vm.count("mvalue_p_thres")) 
		{
			_mvaluePvalueThreshold = vm["mvalue_p_thres"].as<float>();
		}
		if (vm.count("mvalue_method")) 
		{
			_mvalueMethod = vm["mvalue_method"].as<std::string>();
		}
		if (_mvalueMethod == "mcmc") 
		{
			if (vm.count("mcmc_sample")) 
			{
				_mcmcSample = vm["mcmc_sample"].as<unsigned int>();
			}
			if (vm.count("mcmc_burnin")) 
			{
				_mcmcBurnin = vm["mcmc_burnin"].as<unsigned int>();
			}
			if (vm.count("mcmc_prob_random_move")) 
			{
				_mcmcProbRandom = vm["mcmc_prob_random_move"].as<float>();
			}
			if (vm.count("mcmc_max_num_flip")) 
			{
				_mcmcMaxNumFlip = vm["mcmc_max_num_flip"].as<float>();
			}
		}
	}
	if (vm.count("binary_effects")) 
	{
		_willComputeBinaryEffects = true;
		if (vm.count("binary_effects_sample")) 
		{
			_binaryEffectsSample = vm["binary_effects_sample"].as<unsigned int>();
		}
		if (vm.count("binary_effects_large")) 
		{
			_binaryEffectsLargeSample = vm["binary_effects_large"].as<unsigned int>();
		}
		if (vm.count("binary_effects_p_thres")) 
		{
			_binaryEffectsPvalueThreshold = vm["binary_effects_p_thres"].as<float>();
		}
	}
	if (vm.count("seed")) 
	{
		_seed = vm["seed"].as<unsigned int>();
	}
	if (vm["verbose"].as<bool>() == true) 
	{
		_verbose = true;
	}
	if (vm.count("thread")) 
	{
		_threadNum = vm["thread"].as<unsigned int>();
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
	if (_inputFile == "") 
	{
		printErrorAndQuit("A valid input file must be specified using option -input");
	}
	if (_inputLambdaMeanEffects <= 0.0) 
	{
		printErrorAndQuit("lambda_mean option takes a float value > 0");
	}
	if (_inputLambdaHeterogeneity <= 0.0) 
	{
		printErrorAndQuit("lambda_hetero option takes a float value > 0");
	}
	if (_priorSigma <= 0.0) 
	{
		printErrorAndQuit("mvalue_prior_sigma option takes a float value > 0");
	}
	if (_priorAlpha <= 0.0 || _priorBeta <= 0.0) 
	{
		printErrorAndQuit("mvalue_prior_beta option takes two float values > 0");
	}
	if (_mvaluePvalueThreshold < 0.0 || _mvaluePvalueThreshold > 1.0) 
	{
		printErrorAndQuit("mvalue_p_thres takes a float value between 0 and 1");
	}
	if (_mvalueMethod != "exact" &&
		_mvalueMethod != "mcmc" &&
		_mvalueMethod != "variational") 
	{
		printErrorAndQuit("mvalue_method option only takes a value 'exact' or 'mcmc'");
	}
	if (_mcmcSample < 1) 
	{
		printErrorAndQuit("mcmc_sample option takes an integer value > 0");
	}
	if (_mcmcBurnin < 1) 
	{
		printErrorAndQuit("mcmc_burnin option takes an integer value > 0");
	}
	if (_mcmcSample < _mcmcBurnin) 
	{
		printErrorAndQuit("mcmc_sample must be larger than mcmc_burnin");
	}
	if (_mcmcProbRandom < 0.0 || _mcmcProbRandom > 1.0) 
	{
		printErrorAndQuit("mcmc_prob_random takes a float value between 0 and 1");
	}
	if (_mcmcMaxNumFlip <= 0.0) 
	{
		printErrorAndQuit("mcmc_max_num_flip takes a value > 0");
	}
	if (_binaryEffectsSample < 1) 
	{
		printErrorAndQuit("binary_effects_sample option takes an integer value > 0");
	}
	if (_binaryEffectsLargeSample < 1) 
	{
		printErrorAndQuit("binary_effects_large option takes an integer value > 0");
	}
	if (_binaryEffectsPvalueThreshold < 0.0 || _binaryEffectsPvalueThreshold > 1.0) 
	{
		printErrorAndQuit("binary_effects_p_thres takes a float value between 0 and 1");
	}
	
	// Make summary for printing
	
	_argsSummary = "";
	for (int i = 0; i < argc;i++) 
	{
		_argsSummary += std::string(argv[i]) + " ";
	}
	_argsSummary += "\n";
}

void errorHandler(std::string msg)
{
	printf("ERROR: %s\n", msg.c_str());
	std::exit(-1);
}

void doMetaAnalysis()
{
	srand(seed_);
	_numSNPs = _maxNumStudy = 0;
	
	FILE *in, *out;

	try{
		in 	= fopen(_inputFile.c_str(), "r");
		out = fopen(_outputFile.c_str(), "w");

		if(in == NULL || out == NULL)
		{
			throw ERROR_IO_FILE_OPEN_CLOSE;
		}





	}catch(int e){
		printf("Error: [%d] %s\n", e, ERROR_MESSAGE[e]);
	}
	catch(std::exception e){
		printf("Error: unable to Analyze\n");
		exit(ERROR_UNDEFINED);	
	}




	return;
}