#include "MetaSnp.h"
#include "Scheduler.h"
#include <ctime>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <mutex>
#include <condition_variable>

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


// Internally used variables
static std::vector<double> meanEffectParts_;
static std::vector<double> heterogeneityParts_;

const double expectedMedianHanEskinHeterogeneityPart_[] = // from nStudy 2 to 50
{ 0.2195907137,0.2471516439,0.2642270318,0.2780769264,0.2886280267,0.2977812664,0.3020051148,0.3091428179,0.3158605559,0.3221788173,0.3259133140,0.3295976587,0.3335375196,0.3358395088,0.3368309971,0.3421941686,0.3448030927,0.3463590948,0.3477384754,0.3487900288,0.3494938171,0.3542087791,0.3573286353,0.3589703411,0.3586951356,0.3596101209,0.3605611682,0.3624799993,0.3648322669,0.3659817739,0.3671267389,0.3693952373,0.3693395144,0.3696863113,0.3706067524,0.3718103285,0.3749536619,0.3758886239,0.3753612342,0.3781458299,0.3798346038,0.3763434983,0.3796968747,0.3784334922,0.3794411347,0.3808582942,0.3813485882,0.3843230993,0.3824863479 };

//void* thr_func(void* args)
//{
//	thread_struct* thr_args = (thread_struct*)args;
//	std::vector<std::string> LINES = thr_args->lines;
//	FILE* outFile = thr_args->outFile;
//	MetaSnp* metaSnp;    // Store only 1 Snp at a time in memory.
//
//	for (std::string line : LINES)
//	{
//		std::vector<std::string> tokens;
//		split(tokens, line, " ");
//		if (tokens.size() > 1) 
//		{             // only if non-empty
//			if (tokens.at(0).at(0) != '#') 
//			{ // only if non-comment
//				std::string rsid = tokens.at(0);
//				metaSnp = new MetaSnp(rsid);
//				if (tokens.size() % 2 == 0)
//				{
//					printf("WARNING: # of Columns must be odd including Rsid. Last column is ignored.");
//				}
//
//				int nStudy = (int)((tokens.size() - 1) / 2);
//				if (nStudy > maxNumStudy_) 
//				{
//					maxNumStudy_ = nStudy;
//				}
//
//				for (int i = 0; i < nStudy; i++) 
//				{
//					double beta;
//					double standardError;
//					if (tokens.at((2.0 * i) + 1.0).compare("NA") == 0 ||
//						tokens.at((2.0 * i) + 1.0).compare("N/A") == 0 ||
//						tokens.at((2.0 * i) + 1.0).compare("NA") == 0 ||
//						tokens.at((2.0 * i) + 1.0).compare("N/A") == 0) 
//					{
//						metaSnp->addNaStudy();
//					}
//					else 
//					{
//						try {
//							beta = stod(tokens.at((2.0 * i) + 1.0));
//							standardError = stod(tokens.at(2.0 * i + 2.0));
//							if (standardError <= 0.0) 
//							{
//								printf("Standard error cannot be <= zero (%f th column is %f) in the following line.\n",
//									2.0 * i + 3.0, standardError);
//								printf("%s", line.c_str());
//								std::exit(-1);
//							}
//							metaSnp->addStudy(beta, standardError);
//						}
//						catch (std::exception es) 
//						{
//							printf("Incorrect float value in following line. Possibly not a double");
//							printf("%s", line.c_str());
//							std::exit(-1);
//						}
//					}
//				}
//				if (metaSnp->getNStudy() > 1) 
//				{
//					metaSnp->initMem();
//					// Analyze 1 Snp on-the-fly.
//					if (isVerbose_ && numSnps_ % 1000 == 0) 
//					{
//						printf("Analyzing SNP #%d (%s)\n", numSnps_ + 1, rsid.c_str());
//					}
//					// FE, RE, and New RE
//					metaSnp->computeFixedEffects();
//					metaSnp->computeRandomEffects();
//					metaSnp->computeHanEskin(inputLambdaMeanEffect_, inputLambdaHeterogeneity_);
//
//					meanEffectParts_.push_back(metaSnp->getStatisticHanEskinMeanEffectPart());
//					
//					double h = metaSnp->getStatisticHanEskinHeterogeneityPart();
//					
//					if (h > 0.0) 
//					{
//						heterogeneityParts_.push_back(h);
//					}
//					// Binary effects model
//					if (willComputeBinaryEffects_) 
//					{
//						metaSnp->computeBinaryEffectsPvalue(binaryEffectsSample_, rand());
//						if (metaSnp->getPvalueBinaryEffects() <= binaryEffectsPvalueThreshold_) 
//						{
//							metaSnp->computeBinaryEffectsPvalue(binaryEffectsLargeSample_, rand());
//						}
//					}
//					// Mvalues
//					if (willComputeMvalue_) 
//					{
//						if (metaSnp->getPvalueFixedEffects() <= mvaluePvalueThreshold_ ||
//							metaSnp->getPvalueHanEskin() <= mvaluePvalueThreshold_) 
//						{
//							if (mvalueMethod_.compare("exact") == 0) 
//							{
//								metaSnp->computeMvalues(priorAlpha_, priorBeta_, priorSigma_);
//							}
//							else if (mvalueMethod_.compare("mcmc") == 0) 
//							{
//								metaSnp->computeMvaluesMCMC(priorAlpha_, priorBeta_, priorSigma_,
//									mcmcSample_, mcmcBurnin_, mcmcProbRandom_,
//									mcmcMaxNumFlip_,
//									rand());
//							}
//							else 
//							{
//								std::cout << mvalueMethod_ << std::endl;
//								assert(false);
//							}
//						}
//					}
//					numSnps_++;
//				}// end of if
//				mtx.lock();
//				metaSnp->printResults(outFile);
//				mtx.unlock();
//			}// end of if('#')
//		}// end of if(token.size()>1)
//		tokens.clear();
//
//		done_mtx.lock();
//		*thr_args->done_count_ptr = *thr_args->done_count_ptr+1;
//		done_mtx.unlock();
//	}
//}

//void doMetaAnalysis() 
//{
//	srand(seed_);
//	numSnps_ = 0;
//	maxNumStudy_ = 0;
//	meanEffectParts_ = std::vector<double>();
//	heterogeneityParts_ = std::vector<double>();
//	FILE* inFile, *outFile;
//	try {
//		inFile = fopen(inputFile_.c_str(), "r");
//	}
//	catch (std::exception e) 
//	{
//		printErrorAndQuit("ERROR: Input file cannot be opened");
//	}
//	try {
//		outFile = fopen(outputFile_.c_str(), "w");
//	}
//	catch (std::exception e) 
//	{
//		printErrorAndQuit("ERROR: Ouput file cannot be opened");
//	}
//
//	// Print headings
//	MetaSnp::printHeadings(outFile);
//
//	try {
//		std::ifstream inStream(inputFile_);
//		std::vector<std::string> LINES;
//		unsigned int done_count = 0;
//		std::string tmp_ln;
//		thread_struct* args;
//
//#ifdef LINUX
//		int err;
//		pthread_t* thr_ptr = new pthread_t();
//		std::vector<pthread_t*> thread_pool;
//#elif WINDOWS
//		std::thread* thr_ptr;
//		std::vector<std::thread*> thread_pool;
//#endif
//
//		while (std::getline(inStream, tmp_ln))
//		{
//			LINES.push_back(tmp_ln);
//		}
//
//		int N = LINES.size() / threadNum_ + (LINES.size()%threadNum_ != 0);
//
//		//TODO: Divide Lines;
//		for (int nidx = 0; nidx < threadNum_; nidx++)
//		{
//			std::vector<std::string> args_lines;
//			for (int i = N * nidx; i < N * (nidx + 1) && i < LINES.size(); i++) 
//			{
//				args_lines.push_back(LINES.at(i));
//			}
//			args = new thread_struct();
//			args->outFile = outFile;
//			args->done_count_ptr = &done_count;
//			args->lines = args_lines;;
//
//		#ifdef LINUX
//			err = pthread_create(thr_ptr, NULL, &thr_func, args);
//			if (err != 0)
//			{
//				// printf("\nCan't create thread: [%s]", strerr(err));
//				exit(ERROR_THREAD_CREATE);
//			}
//		#elif WINDOWS
//			thr_ptr = new std::thread(thr_func, args);
//		#endif
//
//			thread_pool.push_back(thr_ptr);
//		}
//
//		while (done_count < LINES.size()) {
//			//std::this_thread::sleep_for(std::chrono::microseconds(100));
//			std::cout << "Current Progress : " << done_count <<"/"<<LINES.size()<< " finished." << "\r";
//		}
//
//
//
//
//		// Join Threads
//		for (int nidx = 0; nidx < thread_pool.size(); nidx++)
//		{
//		#ifdef LINUX
//			err = pthread_join(*(pthread_t*)thread_pool.at(nidx), NULL);
//			if (err != 0)
//			{
//				// printf("\nCan't join thread: [%s]\n", strerror(err));
//				exit(ERROR_THREAD_JOIN);
//			}
//		#elif WINDOWS
//			thread_pool.at(nidx)->join();
//		#endif
//		}// end of for(nidx) : Thread Join
//
//
//		try {
//			fclose(inFile);
//			fclose(outFile);
//		}
//		catch (std::exception e)
//		{
//			printf("ERROR: file cannot be closed");
//			std::exit(ERROR_IO_FILE_CLOSE);
//		}
//	}
//	catch (std::exception e) 
//	{
//		printf("ERROR: doMetaAnalysis [%s]", e.what());
//		std::exit(ERROR_META_ANALYSIS);
//	}
//	
//
//	// Reorder the Result File
//	try {
//		std::string tmp_str;
//		FILE* file = fopen(outputFile_.c_str(), "r");
//		std::ifstream Instream(outputFile_.c_str());
//
//		std::getline(Instream, tmp_str); // ignore first line
//		
//		std::vector<ThreadResult> total;
//
//		while (std::getline(Instream, tmp_str)) 
//		{
//			ThreadResult mt(stoi(tmp_str.substr(0, tmp_str.find('\t'))), tmp_str);
//			total.push_back(mt);
//		}
//
//		sort(total.begin(), total.end(), compareKeyValues);
//		Instream.close();
//		FILE* outfile = fopen(outputFile_.c_str(), "w");
//		MetaSnp::printHeadings(outfile);
//		for (int i = 0; i < total.size(); i++) 
//		{
//			fprintf(outfile, "%s\n", total.at(i).values_str.c_str());
//		}
//		
//		fclose(outfile);
//	}
//	catch (std::exception e) 
//	{
//		printf("ERROR: Posterior.txt file has been altered!\n[%s]",e.what());
//		std::exit(-1);
//	}
//}

//void computeLambda() 
//{
//	double median;
//	double expectedMedian;
//	if (meanEffectParts_.size() > 0) 
//	{
//		std::sort(meanEffectParts_.begin(), meanEffectParts_.end());
//		median = meanEffectParts_.at((int)(meanEffectParts_.size() / 2.0));
//		expectedMedian = pow(boost::math::find_location<boost::math::normal>(boost::math::complement(0, 0.25, 1.0)), 2.0);
//		outputLambdaMeanEffect_ = median / expectedMedian;
//	}
//	if (heterogeneityParts_.size() > 0) 
//	{
//		std::sort(heterogeneityParts_.begin(), heterogeneityParts_.end());
//		median = heterogeneityParts_.at((int)(heterogeneityParts_.size() / 2.0));
//		if (maxNumStudy_ > 50) 
//		{
//			expectedMedian = pow(boost::math::find_location<boost::math::normal>(boost::math::complement(0, 0.25, 1.0)), 2.0);
//		}
//		else 
//		{
//			expectedMedian = expectedMedianHanEskinHeterogeneityPart_[maxNumStudy_ - 2];
//		}
//		outputLambdaHeterogeneity_ = median / expectedMedian;
//	}
//}



int main(int argc, char* argv[]) 
{	
	time_t startTime = time(NULL);

	Options op;
	op.handleArguments(argc, argv);
	op.printArguments();

	Scheduler scheduler(&op);
	scheduler.prepare();
	scheduler.run();
	scheduler.cleanUp();
	
	time_t endTime = time(NULL);
	
	op.printLog(difftime(endTime, startTime) / SECONDS_IN_MINUTE);
	
	std::cout<<"----- Finished\n";
	std::cout << "----- Elapsed time: " << difftime(endTime,startTime)/60.0 << " minutes\n";
	return DONE_NORMAL;
}
