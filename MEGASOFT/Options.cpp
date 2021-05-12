#include"Options.h"
#include<condition_variable>
namespace po = boost::program_options;
VOID Options::initialize()
{
	m_outputFile				= DEFAULT_OUTPUT_FILE;
	m_pvalueTableFile			= DEFAULT_PVALUE_TABLE_FILE;
	m_logFile					= DEFAULT_LOG_FILE;
	m_inputLambdaMeanEffect		= DEFAULT_INPUT_LAMBDA_MEAN_EFFECT;
	m_inputLambdaHeterogeneity	= DEFAULT_INPUT_LAMBDA_HETEROGENEITY;
	m_willComputeMvalue			= DEFAULT_WILL_COMPUTE_MVALUE;
	m_priorSigma				= DEFAULT_PRIOR_SIGMA;
	m_priorAlpha				= DEFAULT_PRIOR_ALPHA;
	m_priorBeta					= DEFAULT_PRIOR_BETA;
	m_mvaluePvalueThreshold		= DEFAULT_MVALUE_PVALUE_THRESHOLD;
	m_mvalueMethod				= DEFAULT_MVALUE_METHOD;
	m_mcmcSample				= DEFAULT_MCMC_SAMPLES;
	m_mcmcBurnin				= DEFAULT_MCMC_BURNIN;
	m_mcmcProbRandom			= DEFAULT_MCMC_PROB_RANDOM;
	m_mcmcMaxNumFlip			= DEFAULT_MCMC_MAX_NUM_FLIP;
	m_willComputeBinaryEffects	= DEFAULT_WILL_COMPUTE_BINARY_EFFECTS;
	m_binaryEffectsSample		= DEFAULT_BINARY_EFFECTS_SAMPLES;
	m_binaryEffectsLargeSample	= DEFAULT_BINARY_EFFECTS_LARGE_SAMPLE;
	m_binaryEffectsPvalueThreshold = DEFAULT_BINARY_EFFECTS_PVALUE_THRESHOLD;
	m_seed						= DEFAULT_SEED;
	m_isVerbose					= DEFAULT_VERBOSE;
	m_CPUNumThread				= DEFAULT_CPU_NUM_THREAD;

	m_numSnps		= 0;
	m_maxNumStudy	= 0;
	return;
}
VOID Options::printErrorAndQuit(std::string msg)
{
	printf("%s\n", msg.c_str());
	std::exit(-1);
}
VOID Options::printLog(DOUBLE time)
{
	try {
		FILE* outFile = fopen(m_logFile.c_str(), "w");

		fprintf(outFile, "Arguments: %s ", m_argsSummary.c_str());
		fprintf(outFile, "Input File: %s\n", m_inputFile.c_str());
		fprintf(outFile, "Output File: %s\n", m_outputFile.c_str());
		fprintf(outFile, "Log File: %s\n", m_logFile.c_str());
		fprintf(outFile, "p-value Table File: %s\n", m_pvalueTableFile.c_str());
		fprintf(outFile, "Number of SNPs analyzed: %d\n", m_numSnps);
		fprintf(outFile, "Maximum number of studies: %d\n", m_maxNumStudy);
		fprintf(outFile, "Specified lambda for   mean effect part (default = 1.0): %f\n", m_inputLambdaMeanEffect);
		fprintf(outFile, "Specified lambda for heterogeneity part (default = 1.0): %f\n", m_inputLambdaHeterogeneity);
		fprintf(outFile, "Newly calculated inflation factor lambda for   mean effect part: %f\n", m_outputLambdaMeanEffect);
		fprintf(outFile, "Newly calculated inflation factor lambda for heterogeneity part: %f\n", m_outputLambdaHeterogeneity);
		fprintf(outFile, "Elapsed time : %.3f minutes\n", time);
		fclose(outFile);
	}
	catch (std::exception e)
	{
		printf("ERROR: error encountered while writing in log file");
		std::exit(-1);
	}
}
VOID Options::handleArguments(INT32 argc, char* argv[])
{
	po::options_description desc("Allowed Options");
	desc.add_options()
		("help", "Print help")
		("input", po::value<std::string>(), "Input file (Required)")
		("output", po::value<std::string>(), "Output file (default='out')")
		("pvalue_table", po::value<std::string>(), "value table file (default='HanEskinPvalueTable.txt')")
		("log", po::value<std::string>(), "Log file (default='log')")
		("lambda_mean", po::value<double>(), "(Random Effects) User-specified lambda for mean effect part (default=1.0)")
		("lambda_hetero", po::value<double>(), "(Random Effects) User-specified lambda for heterogeneity part (default=1.0)")
		("mvalue", po::bool_switch()->default_value(false), "Compute m-value(default=FALSE)")
		("mvalue_prior_sigma", po::value<double>(), "Sigma value for normal prior N(0, sigma^2) for effect (default=0.2)")
		("mvalue_prior_alpha", po::value<std::string>(), "Alpha value for Beta dist prior Betadist(alpha,beta) for existence of effect (default=1.0,1.0)")
		("mvalue_prior_beta", po::value<std::string>(), "Beta value for Beta dist prior Betadist(alpha,beta) for existence of effect (default=1.0,1.0)")
		("mvalue_p_thres", po::value<double>(), "Compute m-values only for SNPs whose FE or RE2 p-value is below this threshold (default=1E-7)")
		("mvalue_method", po::value<std::string>(), "Which method to use to calculate m-value between 'exact' and 'mcmc' (default=exact)")
		("mcmc_sample", po::value<long>(), "(MCMC) Number of samples (default=10,000)")
		("mcmc_burnin", po::value<long>(), "(MCMC) Number of burn-in (default=1,000)")
		("mcmc_prob_random_move", po::value<double>(), "(MCMC) Probability that a complete randomization move is suggested (default=0.01)")
		("mcmc_max_num_flip", po::value<double>(), "(MCMC) Usual move is flipping N bits where N ~ U(1,max_num_flip). If an integer value i >= 1 is given, max_num_flip = i. If a float value 0 < k < 1 is given, max_num_flip = k * #studies. (default=0.1)")
		("binary_effects", po::value<std::string>(), "Compute binary effects model p-value (default=FALSE)")
		("binary_effects_sample", po::value<long>(), "(Binary effects) Number of importance sampling samples (default=1,000)")
		("binary_effects_large", po::value<long>(), "(Binary effects) Large number of importance sampling samples for p-values above threshold (default=100,000)")
		("binary_effect_p_thres", po::value<double>(), "(Binary effects) P-value threshold determining if we will use large number of samples (default=1E-4)")
		("seed", po::value<int>(), "Random number generator seed (default=0)")
		("verbose", po::bool_switch()->default_value(false), "Print RSID verbosely per every 1,000 SNPs (default=FALSE)")
		("thread_cpu", po::value<int>(), "Set Number of CPU Threads")
		;

	po::variables_map vm;
	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
	}
	catch (po::error e)
	{
		std::cout << e.what() << std::endl;
		exit(-200);
	}
	if (vm.count("help"))
	{
		std::cout << desc << std::endl;
		std::exit(1);
	}
	if (vm.count("input"))
	{
		m_inputFile = vm["input"].as<std::string>();
	}

	if (vm.count("output"))
	{
		m_outputFile = vm["output"].as<std::string>();
	}
	if (vm.count("pvalue_table"))
	{
		m_pvalueTableFile = vm["pvalue_table"].as<std::string>();
	}
	if (vm.count("log"))
	{
		m_logFile = vm["log"].as<std::string>();
	}
	if (vm.count("lambda_mean"))
	{
		m_inputLambdaMeanEffect = vm["lambda_mean"].as<DOUBLE>();
	}
	if (vm.count("lambda_hetero"))
	{
		m_inputLambdaHeterogeneity = vm["lambda_hetero"].as<DOUBLE>();
	}
	if (vm.count("mvalue"))
	{
		m_willComputeMvalue = TRUE;
		if (vm.count("mvalue_prior_sigma"))
		{
			m_priorSigma = vm["mvalue_prior_sigma"].as<DOUBLE>();
		}
		if (vm.count("mvalue_prior_alpha"))
		{
			m_priorAlpha = stod(vm["mvalue_prior_alpha"].as<std::string>());
		}
		if (vm.count("mvalue_prior_beta"))
		{
			m_priorBeta = stod(vm["mvalue_prior_beta"].as<std::string>());
		}
		if (vm.count("mvalue_p_thres"))
		{
			m_mvaluePvalueThreshold = vm["mvalue_p_thres"].as<DOUBLE>();
		}
		if (vm.count("mvalue_method"))
		{
			m_mvalueMethod = vm["mvalue_method"].as<std::string>();
		}
		if (m_mvalueMethod == "mcmc")
		{
			if (vm.count("mcmc_sample"))
			{
				m_mcmcSample = vm["mcmc_sample"].as<LONG>();
			}
			if (vm.count("mcmc_burnin"))
			{
				m_mcmcBurnin = vm["mcmc_burnin"].as<LONG>();
			}
			if (vm.count("mcmc_prob_random_move"))
			{
				m_mcmcProbRandom = vm["mcmc_prob_random_move"].as<DOUBLE>();
			}
			if (vm.count("mcmc_max_num_flip"))
			{
				m_mcmcMaxNumFlip = vm["mcmc_max_num_flip"].as<DOUBLE>();
			}
		}
	}
	if (vm.count("binary_effects"))
	{
		m_willComputeBinaryEffects = TRUE;
		if (vm.count("binary_effects_sample"))
		{
			m_binaryEffectsSample = vm["binary_effects_sample"].as<LONG>();
		}
		if (vm.count("binary_effects_large"))
		{
			m_binaryEffectsLargeSample = vm["binary_effects_large"].as<LONG>();
		}
		if (vm.count("binary_effects_p_thres"))
		{
			m_binaryEffectsPvalueThreshold = vm["binary_effects_p_thres"].as<DOUBLE>();
		}
	}
	if (vm.count("seed"))
	{
		m_seed = vm["seed"].as<int>();
	}
	if (vm["verbose"].as<bool>() == TRUE)
	{
		m_isVerbose = TRUE;
	}
	if (vm.count("thread_cpu"))
	{
		m_CPUNumThread = vm["thread_cpu"].as<int>();
	}
	if (vm.count("help"))
	{
		std::cout << "------------------------------------------------" << std::endl;
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
	if (m_inputFile== "")
	{
		printErrorAndQuit("A valid input file must be specified using option -input");
	}
	if (m_inputLambdaMeanEffect <= 0.0)
	{
		printErrorAndQuit("lambda_mean option takes a float value > 0");
	}
	if (m_inputLambdaHeterogeneity <= 0.0)
	{
		printErrorAndQuit("lambda_hetero option takes a float value > 0");
	}
	if (m_priorSigma <= 0.0)
	{
		printErrorAndQuit("mvalue_prior_sigma option takes a float value > 0");
	}
	if (m_priorAlpha <= 0.0 || m_priorBeta <= 0.0)
	{
		printErrorAndQuit("mvalue_prior_beta option takes two float values > 0");
	}
	if (m_mvaluePvalueThreshold < 0.0 || m_mvaluePvalueThreshold > 1.0)
	{
		printErrorAndQuit("mvalue_p_thres takes a float value between 0 and 1");
	}
	if (m_mvalueMethod != "exact" &&
		m_mvalueMethod != "mcmc" &&
		m_mvalueMethod != "variational")
	{
		printErrorAndQuit("mvalue_method option only takes a value 'exact' or 'mcmc'");
	}
	if (m_mcmcSample < 1)
	{
		printErrorAndQuit("mcmc_sample option takes an integer value > 0");
	}
	if (m_mcmcBurnin < 1)
	{
		printErrorAndQuit("mcmc_burnin option takes an integer value > 0");
	}
	if (m_mcmcSample < m_mcmcBurnin)
	{
		printErrorAndQuit("mcmc_sample must be larger than mcmc_burnin");
	}
	if (m_mcmcProbRandom < 0.0 || m_mcmcProbRandom > 1.0)
	{
		printErrorAndQuit("mcmc_prob_random takes a float value between 0 and 1");
	}
	if (m_mcmcMaxNumFlip <= 0.0)
	{
		printErrorAndQuit("mcmc_max_num_flip takes a value > 0");
	}
	if (m_binaryEffectsSample < 1)
	{
		printErrorAndQuit("binary_effects_sample option takes an integer value > 0");
	}
	if (m_binaryEffectsLargeSample < 1)
	{
		printErrorAndQuit("binary_effects_large option takes an integer value > 0");
	}
	if (m_binaryEffectsPvalueThreshold < 0.0 || m_binaryEffectsPvalueThreshold > 1.0)
	{
		printErrorAndQuit("binary_effects_p_thres takes a float value between 0 and 1");
	}

	// Make summary for printing
	createSummary(argc, argv);
}
VOID Options::createSummary(INT32 argc, char* argv[])
{
	m_argsSummary = "";
	for (int i = 0; i < argc; i++)
	{
		m_argsSummary += std::string(argv[i]) + " ";
	}
	m_argsSummary += "\n";
}
VOID Options::printArguments(VOID)
{
	std::cout << "Arguments: " + m_argsSummary << std::endl;
}