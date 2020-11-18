#pragma once
#include"MetaSnp.h"
#include"gpu.cuh"
#include"MetaGPU.cuh"
#include <ctime>
#include <boost/program_options.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include<mutex>
#include<condition_variable>
namespace po = boost::program_options;

struct {
	std::vector<std::string> lines;
	FILE* outFile;
	unsigned int* done_count_ptr;
}typedef thread_struct;

// GPU Variable
static int GPU_EXISTS	= false;
static int nBlock		= 0;
static int nThread		= 0;

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
static std::string  mvalueMethod_ = "exact";
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
void handleArguments(int argc, char* argv[]);
void reorderOutputFile();
void* thr_func(void* args);
void runCPU(unsigned int& done_count, std::string& tmp_ln, std::ifstream& inStream, std::vector<std::string>& LINES, thread_struct* args, FILE*& inFile, FILE*& outFile);

bool runGPU(unsigned int& done_count, std::string& tmp_ln, std::ifstream& inStream, std::vector<std::string>& LINES, thread_struct* args, FILE*& inFile, FILE*& outFile);
void doMetaAnalysis();
void computeLambda();
void printLog(double time);
