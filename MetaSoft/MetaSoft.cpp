#include"MetaSnp.h";
#include <ctime>
using namespace std;

// Arguments and default values
static std::string  inputFile_ = "";
static std::string  outputFile_ = "out";
static std::string  logFile_ = "log";
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

void handleArguments(std::string args[]) {
	if (args->size() == 0) {
		cout << "ERROR: No argument. Please type './metasoft -help' to see a list of options" << endl;
		exit(-1);
	}

	for (int i = 0; i < args->size(); i++) {
		cout << args[i] << endl;
	}

	exit(-1);
}

void printErrorAndQuit(std::string msg) {
	printf("%s\n",msg);
	exit(-1);
}

void doMetaAnalysis() {
	srand(time(NULL));
	numSnps_ = 0;
	maxNumStudy_ = 0;
	meanEffectParts_ = std::vector<double>();
	heterogeneityParts_ = std::vector<double>();
	MetaSnp* metaSnp;    // Store only 1 Snp at a time in memory.
	FILE* inFile, *outFile;
	try {
		inFile = fopen(inputFile_.c_str(), "r");
	}
	catch (exception e) {
		printErrorAndQuit("ERROR: Input file cannot be opened");
	}
	try {
		outFile = fopen(outputFile_.c_str(), "w");
	}
	catch (exception e) {
		printErrorAndQuit("ERROR: Ouput file cannot be opened");
	}
	// Print headings
	MetaSnp::printHeadings(outFile);
	
	try {
		std::ifstream inStream(inputFile_);
		// Read 1 Snp information
		std::string readLine;
		while (std::getline(inStream, readLine)) {
			vector<std::string> tokens = split(readLine,"\\s+");
			if (tokens.size() > 1) {             // only if non-empty
				if (tokens[0].at(0) != '#') { // only if non-comment
					std::string rsid = tokens[0];
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
						if (tokens[2 * i + 1].compare("NA") ||
							tokens[2 * i + 1].compare("N/A") ||
							tokens[2 * i + 2].compare("NA") ||
							tokens[2 * i + 2].compare("N/A")) {
							metaSnp->addNaStudy();
						}
						else {
							try {
								beta = stod(tokens[2 * i + 1]);
								standardError = stod(tokens[2 * i + 2]);
								if (standardError <= 0.0) {
									printf("Standard error cannot be <= zero (%d th column is %f) in the following line.\n",
										2 * i + 3, standardError);
									printf("%s",readLine);
									exit(-1);
								}
								metaSnp->addStudy(beta, standardError);
							}
							catch (exception es) {
								printf("Incorrect float value in following line. Possibly not a double");
								printf("%s",readLine);
								exit(-1);
							}
						}
					}
					if (metaSnp->getNStudy() > 1) {
						// Analyze 1 Snp on-the-fly.
						if (isVerbose_ && numSnps_ % 1000 == 0) {
							printf("Analyzing SNP #%d (%s)\n", numSnps_ + 1, rsid);
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
								if (mvalueMethod_.compare("exact")) {
									metaSnp->computeMvalues(priorAlpha_, priorBeta_, priorSigma_);
								}
								else if (mvalueMethod_.compare("mcmc")) {
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
					metaSnp->printResults(outputFile_);

				}
			}
		}
	}
	catch (exception e) {
		printf("ERROR: error encountered while reading input file");
		exit(-1);
	}
	try {
		fclose(inFile);
	}
	catch (exception e) {
		printf("ERROR: file cannot be closed");
		exit(-1);
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
		
		fprintf(outFile,"Arguments: %s ", argsSummary_);
		fprintf(outFile,"Input File: %s\n", inputFile_);
		fprintf(outFile,"Output File: %s\n", outputFile_);
		fprintf(outFile,"Log File: %s\n", logFile_);
		fprintf(outFile,"p-value Table File: %s\n", pvalueTableFile_);
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
		exit(-1);
	}
}


int main(std::string args[]) {
	time_t startTime = time(NULL);
	handleArguments(args);
	printf("Arguments: %s", argsSummary_);
	MetaSnp::readPvalueTableFile(pvalueTableFile_);
	printf("----- Performing meta-analysis\n");
	doMetaAnalysis();
	computeLambda();
	printLog();
	printf("----- Finished\n");
	time_t endTime = time(NULL);
	printf("----- Elapsed time: %.2f minutes\n",
		(endTime - startTime) / (60 * 1000.0F));

	return NORMAL_EXECUTION;
}