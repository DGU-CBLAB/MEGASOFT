#include"MetaSnp.h";
#include <ctime>
#include<boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;

// Multi-Thread variables
static int threadNum_ = 1;
static std::mutex mtx;
// Arguments and default values
static std::string  inputFile_ = "";
static std::string  outputFile_ = "out";
static std::string  logFile_ = ".\log";
static std::string  pvalueTableFile_ = "HanEskinPvalueTable.txt";
static double  inputLambdaMeanEffect_ = 1.0;
static double  inputLambdaHeterogeneity_ = 1.0;
static bool willComputeMvalue_ = false;
static double  priorSigma_ = 0.2;
static double  priorAlpha_ = 1.0;
static double  priorBeta_ = 1.0;
static double  mvaluePvalueThreshold_ = 1E-7;
static std::string  mvalueMethod_ = "exact";
static long    mcmcSample_ = 10000;
static long    mcmcBurnin_ = 1000;
static double  mcmcProbRandom_ = 0.01;
static double  mcmcMaxNumFlip_ = 0.1;
static bool willComputeBinaryEffects_ = false;
static long    binaryEffectsSample_ = 1000;
static long    binaryEffectsLargeSample_ = 100000;
static double  binaryEffectsPvalueThreshold_ = 1E-4;
static int     seed_ = 0;
static bool isVerbose_ = false;
// Internally used variables
static int    numSnps_;
static int    maxNumStudy_;
static double outputLambdaMeanEffect_;
static double outputLambdaHeterogeneity_;
static vector<double> meanEffectParts_;
static vector<double> heterogeneityParts_;
static std::string argsSummary_;

const double expectedMedianHanEskinHeterogeneityPart_[] = // from nStudy 2 to 50
{ 0.2195907137,0.2471516439,0.2642270318,0.2780769264,0.2886280267,0.2977812664,0.3020051148,0.3091428179,0.3158605559,0.3221788173,0.3259133140,0.3295976587,0.3335375196,0.3358395088,0.3368309971,0.3421941686,0.3448030927,0.3463590948,0.3477384754,0.3487900288,0.3494938171,0.3542087791,0.3573286353,0.3589703411,0.3586951356,0.3596101209,0.3605611682,0.3624799993,0.3648322669,0.3659817739,0.3671267389,0.3693952373,0.3693395144,0.3696863113,0.3706067524,0.3718103285,0.3749536619,0.3758886239,0.3753612342,0.3781458299,0.3798346038,0.3763434983,0.3796968747,0.3784334922,0.3794411347,0.3808582942,0.3813485882,0.3843230993,0.3824863479 };

void printErrorAndQuit(std::string msg);

void handleArguments(int argc, char* argv[]) {
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
		("mvalue_prior_beta", po::value<std::string>(), "Alpha and Beta value for Beta dist prior Betadist(alpha,beta) for existence of effect (default=1.0,1.0)")
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
	

	if (vm.count("help")) {
		cout << desc << endl;
		std::exit(1);
	}
	if (vm.count("input")) {
		inputFile_ = vm["input"].as<std::string>();
	}
	
	if (vm.count("output")) {
		outputFile_ = vm["output"].as<std::string>();
	}
	if (vm.count("pvalue_table")) {
		pvalueTableFile_ = vm["pvalue_table"].as<std::string>();
	}
	if (vm.count("log")) {
		logFile_ = vm["log"].as<std::string>();
	}
	if (vm.count("lambda_mean")) {
		inputLambdaMeanEffect_ = vm["lambda_mean"].as<double>();
	}
	if (vm.count("lambda_hetero")) {
		inputLambdaHeterogeneity_ = vm["lambda_hetero"].as<double>();
	}
	if (vm.count("mvalue")) {
		willComputeMvalue_ = true;
		if (vm.count("mvalue_prior_sigma")) {
			priorSigma_ = vm["mvalue_prior_sigma"].as<double>();
		}
		if (vm.count("mvalue_prior_beta")) {
			priorAlpha_ = stod(vm["mvalue_prior_beta"].as<std::string>());
			priorBeta_ = stod(vm["mvalue_prior_beta"].as<std::string>());
		}
		if (vm.count("mvalue_p_thres")) {
			mvaluePvalueThreshold_ = vm["mvalue_p_thres"].as<double>();
		}
		if (vm.count("mvalue_method")) {
			mvalueMethod_ = vm["mvalue_method"].as<std::string>();
		}
		if (mvalueMethod_ == "mcmc") {
			if (vm.count("mcmc_sample")) {
				mcmcSample_ = vm["mcmc_sample"].as<long>();
			}
			if (vm.count("mcmc_burnin")) {
				mcmcBurnin_ = vm["mcmc_burnin"].as<long>();
			}
			if (vm.count("mcmc_prob_random_move")) {
				mcmcProbRandom_ = vm["mcmc_prob_random_move"].as<double>();
			}
			if (vm.count("mcmc_max_num_flip")) {
				mcmcMaxNumFlip_ = vm["mcmc_max_num_flip"].as<double>();
			}
		}
	}
	if (vm.count("binary_effects")) {
		willComputeBinaryEffects_ = true;
		if (vm.count("binary_effects_sample")) {
			binaryEffectsSample_ = vm["binary_effects_sample"].as<long>();
		}
		if (vm.count("binary_effects_large")) {
			binaryEffectsLargeSample_ = vm["binary_effects_large"].as<long>();
		}
		if (vm.count("binary_effects_p_thres")) {
			binaryEffectsPvalueThreshold_ = vm["binary_effects_p_thres"].as<double>();
		}
	}
	if (vm.count("seed")) {
		seed_ = vm["seed"].as<int>();
	}
	if (vm["verbose"].as<bool>() == true) {
		isVerbose_ = true;
	}
	if (vm.count("thread")) {
		threadNum_ = vm["thread"].as<int>();
	}
	if (vm.count("help")) {
		std::cout << "------------------------------------------------"<<endl;
		std::cout << "The format of input_file:" << endl;
		std::cout << "  Each row is each SNP." << endl;
		std::cout << "  1st column is RSID." << endl;
		std::cout << "  2nd and 3rd column are beta and its standard error for 1st study." << endl;
		std::cout << "  4th and 5th column are beta and its standard error for 2nd study." << endl;
		std::cout << "  6th and 7th column are beta and its standard error for 3rd study." << endl;
		std::cout << "  and so on..." << endl;
		std::cout << "------------------------------------------------" << endl;
		std::cout << endl;
		std::exit(-1);
	}


	std::cout << "-------- Processing arguments ---------" << std::endl;
	if (inputFile_ == "") {
		printErrorAndQuit("A valid input file must be specified using option -input");
	}
	if (inputLambdaMeanEffect_ <= 0.0) {
		printErrorAndQuit("lambda_mean option takes a float value > 0");
	}
	if (inputLambdaHeterogeneity_ <= 0.0) {
		printErrorAndQuit("lambda_hetero option takes a float value > 0");
	}
	if (priorSigma_ <= 0.0) {
		printErrorAndQuit("mvalue_prior_sigma option takes a float value > 0");
	}
	if (priorAlpha_ <= 0.0 || priorBeta_ <= 0.0) {
		printErrorAndQuit("mvalue_prior_beta option takes two float values > 0");
	}
	if (mvaluePvalueThreshold_ < 0.0 || mvaluePvalueThreshold_ > 1.0) {
		printErrorAndQuit("mvalue_p_thres takes a float value between 0 and 1");
	}
	if (mvalueMethod_ != "exact" &&
		mvalueMethod_ != "mcmc" &&
		mvalueMethod_ != "variational") {
		printErrorAndQuit("mvalue_method option only takes a value 'exact' or 'mcmc'");
	}
	if (mcmcSample_ < 1) {
		printErrorAndQuit("mcmc_sample option takes an integer value > 0");
	}
	if (mcmcBurnin_ < 1) {
		printErrorAndQuit("mcmc_burnin option takes an integer value > 0");
	}
	if (mcmcSample_ < mcmcBurnin_) {
		printErrorAndQuit("mcmc_sample must be larger than mcmc_burnin");
	}
	if (mcmcProbRandom_ < 0.0 || mcmcProbRandom_ > 1.0) {
		printErrorAndQuit("mcmc_prob_random takes a float value between 0 and 1");
	}
	if (mcmcMaxNumFlip_ <= 0.0) {
		printErrorAndQuit("mcmc_max_num_flip takes a value > 0");
	}
	if (binaryEffectsSample_ < 1) {
		printErrorAndQuit("binary_effects_sample option takes an integer value > 0");
	}
	if (binaryEffectsLargeSample_ < 1) {
		printErrorAndQuit("binary_effects_large option takes an integer value > 0");
	}
	if (binaryEffectsPvalueThreshold_ < 0.0 || binaryEffectsPvalueThreshold_ > 1.0) {
		printErrorAndQuit("binary_effects_p_thres takes a float value between 0 and 1");
	}
	// Make summary for printing
	
	argsSummary_ = "";
	for (int i = 0; i < argc;i++) {
		argsSummary_ += std::string(argv[i]) + " ";
	}
	argsSummary_ += "\n";
}

void printErrorAndQuit(std::string msg) {
	printf("%s\n",msg.c_str());
	std::exit(-1);
}
void thr_func(std::string readLine, FILE* outFile) {
	MetaSnp* metaSnp;    // Store only 1 Snp at a time in memory.
	vector<std::string> tokens;
	split(tokens, readLine, " ");
	if (tokens.size() > 1) {             // only if non-empty
		if (tokens.at(0).at(0) != '#') { // only if non-comment
			std::string rsid = tokens.at(0);
			metaSnp = new MetaSnp(rsid);
			if (tokens.size() % 2 == 0)
				printf("WARNING: # of Columns must be odd including Rsid. Last column is ignored.");
			int nStudy = (int)((tokens.size() - 1) / 2);
			if (nStudy > maxNumStudy_) {
				maxNumStudy_ = nStudy;
			}
			for (int i = 0; i < nStudy; i++) {
				double beta;
				double standardError;
				if (tokens.at(2 * i + 1).compare("NA") == 0 ||
					tokens.at(2 * i + 1).compare("N/A") == 0 ||
					tokens.at(2 * i + 2).compare("NA") == 0 ||
					tokens.at(2 * i + 2).compare("N/A") == 0) {
					metaSnp->addNaStudy();
				}
				else {
					try {
						beta = stod(tokens.at(2 * i + 1));
						standardError = stod(tokens.at(2 * i + 2));
						if (standardError <= 0.0) {
							printf("Standard error cannot be <= zero (%d th column is %f) in the following line.\n",
								2 * i + 3, standardError);
							printf("%s", readLine.c_str());
							std::exit(-1);
						}
						metaSnp->addStudy(beta, standardError);
					}
					catch (exception es) {
						printf("Incorrect float value in following line. Possibly not a double");
						printf("%s", readLine.c_str());
						std::exit(-1);
					}
				}
			}
			if (metaSnp->getNStudy() > 1) {
				// Analyze 1 Snp on-the-fly.
				if (isVerbose_ && numSnps_ % 1000 == 0) {
					printf("Analyzing SNP #%d (%s)\n", numSnps_ + 1, rsid.c_str());
				}
				// FE, RE, and New RE
				metaSnp->computeFixedEffects();
				metaSnp->computeRandomEffects();
				metaSnp->computeHanEskin(inputLambdaMeanEffect_,
					inputLambdaHeterogeneity_);
				meanEffectParts_.push_back(metaSnp->getStatisticHanEskinMeanEffectPart());
				double h = metaSnp->getStatisticHanEskinHeterogeneityPart();
				if (h > 0.0) {
					heterogeneityParts_.push_back(h);
				}
				// Binary effects model
				if (willComputeBinaryEffects_) {
					metaSnp->computeBinaryEffectsPvalue(binaryEffectsSample_, rand());
					if (metaSnp->getPvalueBinaryEffects() <= binaryEffectsPvalueThreshold_) {
						metaSnp->computeBinaryEffectsPvalue(binaryEffectsLargeSample_, rand());
					}
				}
				// Mvalues
				if (willComputeMvalue_) {
					if (metaSnp->getPvalueFixedEffects() <= mvaluePvalueThreshold_ ||
						metaSnp->getPvalueHanEskin() <= mvaluePvalueThreshold_) {
						if (mvalueMethod_.compare("exact") == 0) {
							metaSnp->computeMvalues(priorAlpha_, priorBeta_, priorSigma_);
						}
						else if (mvalueMethod_.compare("mcmc") == 0) {
							metaSnp->computeMvaluesMCMC(priorAlpha_, priorBeta_, priorSigma_,
								mcmcSample_, mcmcBurnin_, mcmcProbRandom_,
								mcmcMaxNumFlip_,
								rand());
						}
						else {
							std::cout << mvalueMethod_ << endl;
							assert(false);
						}
					}
				}
				numSnps_++;
			}
			mtx.lock();
			metaSnp->printResults(outFile);
			mtx.unlock();
		}
	}
	tokens.clear();
}
void doMetaAnalysis() {
	srand(time(NULL));
	numSnps_ = 0;
	maxNumStudy_ = 0;
	meanEffectParts_ = std::vector<double>();
	heterogeneityParts_ = std::vector<double>();
	FILE* inFile, *outFile;
	try {
		inFile = fopen(inputFile_.c_str(), "r");
	}
	catch (exception e) {
		printErrorAndQuit("ERROR: Input file cannot be opened");
		exit(-1);
	}
	try {
		outFile = fopen(outputFile_.c_str(), "w");
	}
	catch (exception e) {
		printErrorAndQuit("ERROR: Ouput file cannot be opened");
		exit(-1);
	}
	// Print headings
	MetaSnp::printHeadings(outFile);

	// Thread
	try {
		std::ifstream inStream(inputFile_);
		std::string readLine;
		int count = 0;
		std::vector<std::thread> tr_vec;
		while (std::getline(inStream, readLine)) {
			tr_vec.push_back(std::thread(thr_func, readLine, outFile));
			bool b = false;
			while (true) {
				if (b == true || tr_vec.size() < threadNum_) {
					break;
				}
				else {
					std::this_thread::sleep_for(std::chrono::seconds(1));
				}
				for (int k = 0; k < tr_vec.size(); k++) {
					if (tr_vec.at(k).joinable()) {
						tr_vec.at(k).join();
						tr_vec.erase(tr_vec.begin() + k);
						b = true;
						cout << "Current Progress : " << ++count << " finished." << "\r";
						break;
					}
				}
			}
		}
		for (int i = 0; i < tr_vec.size(); i++) {
			tr_vec.at(i).join();
		}
	}
	catch (exception e) {
		printf("ERROR: error encountered while reading input file");
		std::exit(-1);
	}
	try {
		fclose(inFile);
		fclose(outFile);
	}
	catch (exception e) {
		printf("ERROR: file cannot be closed");
		std::exit(-1);
	}

	// Reorder the Result File
	try {
		std::string readLine;
		FILE* file = fopen(outputFile_.c_str(), "r");
		std::ifstream Instream(file);

		std::getline(Instream, readLine); // ignore first line
		
		std::vector<map_tuple> total;

		while (std::getline(Instream, readLine)) {
			cout << readLine.substr(0, readLine.find('\t')) << endl;
			map_tuple mt(stoi(readLine.substr(0, readLine.find('\t'))), readLine);
			total.push_back(mt);
		}

		sort(total.begin(), total.end(), map_comp);
		Instream.close();
		FILE* outfile = fopen(outputFile_.c_str(), "w");
		MetaSnp::printHeadings(outfile);
		cout << total.size() << endl;
		for (int i = 0; i < total.size(); i++) {
			fprintf(outfile, "%s\n", total.at(i).val.c_str());
		}
		
		fclose(outfile);
	}
	catch (exception e) {
		printf("ERROR: Posterior.txt file has been altered!");
		std::exit(-1);
	}
}

void computeLambda() {
	double median;
	double expectedMedian;
	if (meanEffectParts_.size() > 0) {
		std::sort(meanEffectParts_.begin(), meanEffectParts_.begin() + meanEffectParts_.size());
		median = meanEffectParts_.at((int)(meanEffectParts_.size() / 2.0));
		expectedMedian = pow(boost::math::find_location<boost::math::normal>(boost::math::complement(0, 0.25, 1.0)), 2.0);
		outputLambdaMeanEffect_ = median / expectedMedian;
	}
	if (heterogeneityParts_.size() > 0) {
		std::sort(heterogeneityParts_.begin(), heterogeneityParts_.begin() + heterogeneityParts_.size());
		median = heterogeneityParts_.at((int)(heterogeneityParts_.size() / 2.0));
		if (maxNumStudy_ > 50) {
			expectedMedian = pow(boost::math::find_location<boost::math::normal>(boost::math::complement(0, 0.25, 1.0)), 2.0);
		}
		else {
			expectedMedian = expectedMedianHanEskinHeterogeneityPart_[maxNumStudy_ - 2];
		}
		outputLambdaHeterogeneity_ = median / expectedMedian;
	}
}

void printLog() {
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
		fclose(outFile);
	}
	catch (exception e) {
		printf("ERROR: error encountered while writing in log file");
		std::exit(-1);
	}
}


int main(int argc, char* argv[]) {	
	time_t startTime = time(NULL);
	handleArguments(argc, argv);
	std::cout << "Arguments: " + argsSummary_ << std::endl;
	MetaSnp::readPvalueTableFile(pvalueTableFile_);
	std::cout<<"----- Performing meta-analysis\n";
	doMetaAnalysis();
	cout << "---- Performing lambda compute\n";
	computeLambda();
	cout << "---- print Log\n";
	printLog();
	std::cout<<"----- Finished\n";
	time_t endTime = time(NULL);
	std::cout << "----- Elapsed time: " << difftime(endTime,startTime)/60.0 << " minutes\n";

	return NORMAL_EXECUTION;
}