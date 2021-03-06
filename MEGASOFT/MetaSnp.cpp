#include"MetaSnp.h"
#include<boost/math/special_functions/beta.hpp>
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
/*
	Functional Replacement
	java : Probability.chiSquareComplemented(v, x)
	c++ : 1 - boost::math::gamma_p<double,double>(v/2.0, x/2.0)
*/
/*
	https://stackoverflow.com/questions/13042561/fatal-error-lnk1104-cannot-open-file-libboost-system-vc110-mt-gd-1-51-lib
	boost LNK problem solved
*/
#define POSINF INFINITY
#define NEGINF -INFINITY

double MetaSnp::ML_ESTIMATE_CHANGE_RATIO_THRESHOLD = 0.00001;
double MetaSnp::LOG_SQRT2PI = 0.5 * log(2 * (double)M_PI);

double MetaSnp::TABLE_MAX_THRESHOLD = 33.0;
double** MetaSnp::pvalueTable_;
bool MetaSnp::isPvalueTableRead_ = false;



double MetaSnp::logBeta(double m, double n) 
{
	return log(boost::math::beta(m,n));
}
double MetaSnp::chiSquareComplemented(double v, double x) 
{
	return (1 - boost::math::gamma_p<double, double>(v / 2.0, x / 2.0));
}
/// <summary>
/// Compare ThreadResult Key Value(for sort algorithm)
/// </summary>
/// <param name="a">Object a</param>
/// <param name="b">Object b</param>
/// <returns>
/// Return 'true' if b.key value is bigger than a.key.
/// Return 'false' if a.key value is bigger than b.key.
/// </returns>
/// <example>
/// <code>
/// sort(vec.begin(), vec.end(), compareKeyValues);
/// </code>
/// </example>
bool compareKeyValues(const ThreadResult& a, const ThreadResult& b) 
{
	return a.key < b.key;
}

/// <summary>
/// Split String into tokens using specified delimiter.
/// </summary>
/// <param name="tokens">String vector to save Results</param>
/// <param name="str">Target String to be splited</param>
/// <param name="delim">delimiter</param>
void split(std::vector<std::string>& tokens, const std::string& str, const std::string& delim)
{
	size_t prev = 0, pos = 0;
	// Empty token results
	tokens.clear();

	do
	{
		pos = str.find(delim, prev);

		if (str.find("\t", prev) < pos) 
		{
			pos = str.find("\t", prev);
		}
		
		if (pos == std::string::npos)
		{
			pos = str.length();
		}

		std::string token = str.substr(prev, pos - prev);
		
		if (!token.empty()) 
		{
			tokens.push_back(token);
		}

		prev = pos + delim.length();

	} while (pos < str.length() && prev < str.length());

	return;
}

/// <summary>
/// MetaSnp Class Constructor
/// </summary>
/// <param name="rsid">Rapid Stain Identification Series (RSID)</param>
MetaSnp::MetaSnp(std::string rsid) 
{
	rsid_			= rsid;
	betas_			= std::vector<double>();
	standardErrors_ = std::vector<double>();
	isNa_			= std::vector<bool>();
	hvalues_		= std::vector<double>();
	mvalues_		= std::vector<double>();
}

/// <summary>
/// Add Study to SNP
/// </summary>
/// <param name="beta">Beta value</param>
/// <param name="standardError">Standard Error</param>
void MetaSnp::addStudy(double beta, double standardError) 
{
	nStudy_++;
	betas_.push_back(beta);
	standardErrors_.push_back(standardError);
	nStudyIncludingNa_++;
	isNa_.push_back(false);
}

/// <summary>
/// Add NA Study to SNP (no parameters required)
/// </summary>
void MetaSnp::addNaStudy() 
{
	nStudyIncludingNa_++;
	isNa_.push_back(true);
}

/// <summary>
/// Compute Fixed Effects
/// </summary>
/// <param name="lambdaMeanEffect">Lambda Mean Effect value</param>
/// <remarks>
/// This Function require memory for betas/variances/weights from MetaSnp::initMen()
/// </remarks>
void MetaSnp::computeFixedEffects(double lambdaMeanEffect) 
{
	for (int i = 0; i < nStudy_; i++) 
	{
		betas[i]		= betas_.at(i);
		variances[i]	= pow(standardErrors_.at(i), 2.0);
		weights[i]		= 1.0 / variances[i];
	}

	double sumBetaProdWeight	= 0.0;
	double sumWeight			= 0.0;

	for (int i = 0; i < nStudy_; i++) 
	{
		sumBetaProdWeight	+= betas[i] * weights[i];
		sumWeight			+= weights[i];
	}

	betaFixedEffects_			= sumBetaProdWeight / sumWeight;
	standardErrorFixedEffects_	= 1.0				/ sqrt(sumWeight);
	statisticFixedEffects_		= sumBetaProdWeight / sqrt(sumWeight);
	statisticFixedEffects_		/= sqrt(lambdaMeanEffect);
	pvalueFixedEffects_			= chiSquareComplemented(1.0, pow(statisticFixedEffects_,2.0));
	isFixedEffectsComputed_		= true;
	
	return;
}

/// <summary>
/// Compute Fixed Effects(default Lambda Mean Effect = 1.0)
/// </summary>
/// <remarks>
/// This function use MetaSnp::computeFixedEffects(double lambdaMeanEffect).
/// </remarks>
void MetaSnp::computeFixedEffects() 
{
	computeFixedEffects(1.0);
	return;
}

/// <summary>
/// Compute Heterogeneity
/// </summary>
/// <remarks>
/// This Function require memory for betas/variances/weights from MetaSnp::initMen()
/// </remarks>
void MetaSnp::computeHeterogeneity() 
{

	for (int i = 0; i < nStudy_; i++) 
	{
		betas[i]		= betas_.at(i);
		variances[i]	= pow(standardErrors_.at(i), 2.0);
		weights[i]		= 1.0 / variances[i];
	}

	double sumBetaProdWeight	= 0.0;
	double sumWeight			= 0.0;
	double sumWeightSquare		= 0.0;
	
	for (int i = 0; i < nStudy_; i++) 
	{
		sumBetaProdWeight	+= betas[i] * weights[i];
		sumWeight			+= weights[i];
		sumWeightSquare		+= weights[i] * weights[i];
	}

	// Compute Q and ISquare
	double meanBeta = sumBetaProdWeight / sumWeight;
	double Q = 0.0;
	
	for (int i = 0; i < nStudy_; i++) 
	{
		Q += weights[i] * (betas[i] - meanBeta) * (betas[i] - meanBeta);
	}
	
	statisticQ_ = Q;
	pvalueQ_	= chiSquareComplemented(nStudy_ - 1.0, Q);
	
	// replace java Math.max : finding maximum value of two
	statisticISquare_ = (Q - (double)(nStudy_ - 1.0))/Q*100.0;
	
	if (statisticISquare_ < 0.0) 
	{
		statisticISquare_ = 0.0;
	}

	// Compute tauSquare
	double meanWeight = sumWeight / nStudy_;
	double Sw2 = (1.0 / (double)(nStudy_ - 1.0))
				* (sumWeightSquare - nStudy_ * meanWeight * meanWeight);
	double U = (double)(nStudy_ - 1.0)
				* (meanWeight - Sw2 / (nStudy_ * meanWeight));
	statisticTauSquare_ = (Q - (double)(nStudy_ - 1.0)) / U;
	
	if (statisticTauSquare_ < 0.0) 
	{
		statisticTauSquare_ = 0.0;
	}
	
	isHeterogeneityComputed_ = true;
	return;
}

/// <summary>
/// Compute Random Effects
/// </summary>
/// <remarks>
/// This Function require memory for betas/variances/weights from MetaSnp::initMen()
/// </remarks>
void MetaSnp::computeRandomEffects() 
{
	if (!isHeterogeneityComputed_)
	{
		computeHeterogeneity();
	}

	for (int i = 0; i < nStudy_; i++) 
	{
		betas[i]		= betas_.at(i);
		variances[i]	= pow(standardErrors_.at(i), 2.0);
		weights[i]		= 1.0 / variances[i];
	}

	double sumBetaProbWeightWithTau = 0.0;
	double sumWeightWithTau = 0.0;
	
	for (int i = 0; i < nStudy_; i++) 
	{
		sumBetaProbWeightWithTau	+= betas[i] / (variances[i] + statisticTauSquare_);
		sumWeightWithTau			+=		1.0 / (variances[i] + statisticTauSquare_);
	}

	betaRandomEffects_			= sumBetaProbWeightWithTau	/ sumWeightWithTau;
	standardErrorRandomEffects_ = 1.0						/ sqrt(sumWeightWithTau);
	statisticRandomEffects_		= sumBetaProbWeightWithTau	/ sqrt(sumWeightWithTau);
	pvalueRandomEffects_		= chiSquareComplemented(1.0, pow(statisticRandomEffects_, 2.0));
	isRandomEffectsComputed_	= true;

	return;
}

/// <summary>
/// Compute Mvalues (Pure Computation)
/// </summary>
/// <param name="priorAlpha">Set Prior Alpha for Beta Distribution</param>
/// <param name="priorBeta">Set Prior Beta for Beta Distribution</param>
/// <param name="priorSigma">Set Prior Sigma</param>
void MetaSnp::computeMvalues(double priorAlpha, double priorBeta, double priorSigma) 
{
	mvalues_.clear();
	double priorVar = priorSigma * priorSigma;
	double* ts		= (double*)malloc(sizeof(double)*nStudy_); // Precision

	for (int i = 0; i < nStudy_; i++) 
	{
		betas[i]	= betas_.at(i);
		ts[i]		= 1.0 / pow(standardErrors_.at(i), 2.0);
	}

	bool* H1			= (bool*)malloc(sizeof(bool)*nStudy_);	// Array defining configuration
	double* priorConfig = (double*)malloc(sizeof(double)*(nStudy_ + 1));	// Prob of each configuration with respect to # of studies with effect

	for (int i = 0; i <= nStudy_; i++) 
	{
		priorConfig[i] = exp(logBeta(i + priorAlpha, nStudy_ - i + priorBeta))
					/ exp(logBeta(priorAlpha, priorBeta));
	}

	double* accumProbH0 = (double*)malloc(sizeof(double)*nStudy_);	// Accumulate probability for each study
	double* accumProbH1 = (double*)malloc(sizeof(double)*nStudy_);

	for (int i = 0; i < nStudy_; i++) 
	{
		accumProbH0[i] = 0.0;
		accumProbH1[i] = 0.0;
	}

	for (int T = 0; T < pow(2, nStudy_); T++) 
	{	// For all configurations.
		int t = T;
		int numH1 = 0;
		for (int i = 0; i < nStudy_; i++) 
		{
			if (t % 2 == 1) 
			{
				H1[i] = true;
				numH1++;
			}
			else 
			{
				H1[i] = false;
			}

			t = (int)(floor(t / 2.0));
		}

		// First for null points
		double probNullPoints = 1.0;
		for (int i = 0; i < nStudy_; i++) 
		{
			if (!H1[i]) 
			{
				probNullPoints *= sqrt(ts[i] / (2 * M_PI)) * exp(-ts[i] * betas[i] * betas[i] / 2);
			}
		}

		// Second for alternative points
		double probAltPoints = 1.0;
		if (numH1 > 0) 
		{
			double sum_t	= 0.0;
			double sum_tm	= 0.0;
			double sum_tmm	= 0.0;
			double prod_t	= 1.0;

			for (int i = 0; i < nStudy_; i++) 
			{
				if (H1[i]) 
				{
					sum_t	+= ts[i];
					sum_tm	+= ts[i] * betas[i];
					sum_tmm += ts[i] * betas[i] * betas[i];
					prod_t	*= ts[i];
				}
			}

			double betaJoint	= sum_tm / sum_t;
			double tJoint		= sum_t;                                
			double tconst		= 1 / ((1 / tJoint) + priorVar);
			double scaleFactor	= sqrt(prod_t / sum_t) * pow(2 * M_PI, -0.5*(numH1 - 1))
								* exp(-(sum_tmm - sum_tm * sum_tm / sum_t) / 2);
			double jointPDF		= sqrt(tconst / (2 * M_PI)) * exp(-tconst * betaJoint * betaJoint / 2);
			probAltPoints		= jointPDF * scaleFactor;
		} // end of if

		for (int i = 0; i < nStudy_; i++) 
		{
			if (H1[i]) 
			{
				accumProbH1[i] += probNullPoints * probAltPoints * priorConfig[numH1];
			}
			else 
			{
				accumProbH0[i] += probNullPoints * probAltPoints * priorConfig[numH1];
			}
		}
	} // end of for(T)

	for (int i = 0; i < nStudy_; i++) 
	{
		double mvalue = accumProbH1[i] / (accumProbH0[i] + accumProbH1[i]);
		mvalues_.push_back(mvalue);
	}

	isMvaluesComputed_ = true;

	free(ts);
	free(H1);
	free(priorConfig);
	free(accumProbH0);
	free(accumProbH1);

	return;
}

/// <summary>
/// Compute Mvalues (Using MCMC)
/// </summary>
/// <param name="priorAlpha">Set Prior Alpha for Beta Distribution</param>
/// <param name="priorBeta">Set Prior Beta for Beta Distribution</param>
/// <param name="priorSigma">Set Prior Sigma</param>
/// <param name="sample">Set Number of Samples for MCMC</param>
/// <param name="burnin">Set burninCount for MCMC</param>
/// <param name="probRandom">Set Probability Random Value MCMC</param>
/// <param name="maxNumFlipArg">Set Maximum Number of Flips for MCMC</param>
/// <param name="seed">Set Seed Number for Random function</param>
void MetaSnp::computeMvaluesMCMC(double priorAlpha, double priorBeta, double priorSigma, long sample, long burnin, double probRandom, double maxNumFlipArg, unsigned int seed) 
{
	srand(seed);

	int maxNumFlip = 1;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> randDouble(0.0, 1.0);

	if (maxNumFlipArg < 1.0) 
	{
		maxNumFlip = (int)floor(maxNumFlipArg * nStudy_);
		if (maxNumFlip < 1)
		{
			maxNumFlip = 1;
		}
	}
	else maxNumFlip = (int)floor(maxNumFlipArg);
	
	mvalues_.clear();
	double priorVar = priorSigma * priorSigma;
	
	double* ts		= (double*)malloc(sizeof(double)*nStudy_); // Precision

	for (int i = 0; i < nStudy_; i++) 
	{
		betas[i] = betas_.at(i);
		ts[i] = 1.0 / pow(standardErrors_.at(i), 2.0);
	}
	
	bool* H1				= (bool*)malloc(sizeof(bool)*nStudy_); // Array defining configuration
	double* logPriorConfig	= (double*)malloc(sizeof(double)*(nStudy_ + 1.0)); // Prob of each configuration with respect to # of studies with effect

	for (int i = 0; i <= nStudy_; i++) 
	{
		logPriorConfig[i] = logBeta(i + priorAlpha, nStudy_ - (double)i + priorBeta)
							- logBeta(priorAlpha, priorBeta);
	}

	int* accumCntH0 = (int*)malloc(sizeof(int)*nStudy_); // Accumullate count for each study
	int* accumCntH1 = (int*)malloc(sizeof(int)*nStudy_);

	for (int i = 0; i < nStudy_; i++) 
	{
		accumCntH0[i] = 0;
		accumCntH1[i] = 0;
	}

	// Start from a random configuration
	int numH1 = 0;
	for (int i = 0; i < nStudy_; i++) 
	{
		H1[i] = (0 + (rand() % (1 - 0 + 1)) == 1);
		if (H1[i] == 1) numH1++;
	}
	
	long burninCount	= burnin;
	long chainCount		= 0;
	bool* tmp			= (bool*)malloc(sizeof(bool)*nStudy_);
	int* shuffleBuffer	= (int*)malloc(sizeof(int)*nStudy_);

	for (int i = 0; i < nStudy_; i++) 
	{
		shuffleBuffer[i] = i;
	}

	// Chain
	while (chainCount < sample) 
	{
		double currentLogProb = observationLogLikelihood(betas, nStudy_, ts, H1, numH1, priorVar) + logPriorConfig[numH1];
		
		if( randDouble(gen) > probRandom)
		{
			// Usual jump
			int numFlip = (rand() % maxNumFlip) + 1;
			if (numFlip > nStudy_) 
			{
				numFlip = numFlip % nStudy_;
			}

			for (int i = 0; i < numFlip; i++) 
			{
				int pick = rand() % (nStudy_ - i);
				std::swap(shuffleBuffer[i], shuffleBuffer[i + pick]);
			}

			for (int i = 0; i < numFlip; i++) 
			{
				int j = shuffleBuffer[i];
				H1[j] = !H1[j];
				numH1 += H1[j] ? 1 : -1;
			}

			double nextLogProb = observationLogLikelihood(betas, nStudy_, ts, H1, numH1, priorVar) + logPriorConfig[numH1];

			if (nextLogProb > currentLogProb || randDouble(gen) < exp(nextLogProb - currentLogProb)) 
			{
				currentLogProb = nextLogProb;
			}
			else 
			{
				// Stay ... revert back
				for (int i = 0; i < numFlip; i++) 
				{
					int j = shuffleBuffer[i];
					H1[j] = !H1[j];
					numH1 += H1[j] ? 1 : -1;
				}
			}
		}
		else 
		{
			// Randomization move
			int tmpNumH1 = 0;
			for (int i = 0; i < nStudy_; i++) 
			{
				tmp[i] = (0 + (rand() % (1 - 0 + 1)) == 1);//rand() % 2;
				if (tmp[i]) 
				{
					tmpNumH1++;
				}
			}

			double nextLogProb = observationLogLikelihood(betas, nStudy_, ts, tmp, tmpNumH1, priorVar) + logPriorConfig[tmpNumH1];
			
			if (nextLogProb > currentLogProb || randDouble(gen) < exp(nextLogProb - currentLogProb)) 
			{
				// Move
				// Copy tmp to H1
				for (int cpy = 0; cpy < nStudy_; cpy++) 
				{
					H1[cpy] = tmp[cpy];
				}

				numH1 = tmpNumH1;
				currentLogProb = nextLogProb;
			
			}
			else 
			{
				// Stay...
			}
		}

		// Are we still in Burn-in?
		if (burninCount > 0) 
		{
			burninCount--;
		}
		else {
		
			for (int i = 0; i < nStudy_; i++) 
			{
				if (H1[i] == true) 
				{
					accumCntH1[i]++;
				}
				else 
				{
					accumCntH0[i]++;
				}
			}
			chainCount++;
		}// end of if
	}// end of while
	
	for (int i = 0; i < nStudy_; i++) 
	{
		double mvalue = (double)accumCntH1[i] / (accumCntH0[i] + accumCntH1[i]);
		mvalues_.push_back(mvalue);
	}
	isMvaluesComputed_ = true;

	free(ts);
	free(H1);
	free(logPriorConfig);
	free(accumCntH0);
	free(accumCntH1);
	free(tmp);
	free(shuffleBuffer);

	return;
}

/// <summary>
/// Calculate Observation Log Likelihood
/// </summary>
/// <param name="betas">Betas for log Probability</param>
/// <param name="betas_size">Size of Betas</param>
/// <param name="ts">???????</param>
/// <param name="H1">For Randomized Movement</param>
/// <param name="numH1">Number of H1</param>
/// <param name="priorVar">Priority variable</param>
/// <returns></returns>
double MetaSnp::observationLogLikelihood(double* betas, int betas_size, double* ts, bool* H1, int numH1, double priorVar) 
{
	int n = betas_size;

	// First for null points
	double logProbNullPoints = 0;
	for (int i = 0; i < n; i++) 
	{
		if (H1[i] == false) 
		{
			logProbNullPoints += 0.5 * log(ts[i]) - LOG_SQRT2PI - ts[i] * betas[i] * betas[i] / 2;
		}
	}// end of for

	// Second for alternative points
	double logProbAltPoints = 0;
	
	if (numH1 > 0) 
	{
		double sum_t	= 0.0;
		double sum_tm	= 0.0;
		double sum_tmm	= 0.0;
		double sum_logt = 0.0;

		for (int i = 0; i < n; i++) 
		{
			if (H1[i]) 
			{
				sum_t		+= ts[i];
				sum_tm		+= ts[i] * betas[i];
				sum_tmm		+= ts[i] * betas[i] * betas[i];
				sum_logt	+= log(ts[i]);
			}
		}

		double betaJoint		= sum_tm / sum_t;
		double tJoint			= sum_t;
		double tconst			= 1 / ((1 / tJoint) + priorVar);
		double logScaleFactor	= 
			-(numH1 - 1)*LOG_SQRT2PI + 0.5*sum_logt - 0.5* log(sum_t)
			- (sum_tmm - sum_tm * sum_tm / sum_t) / 2;
		double logJointPDF		= 0.5 * log(tconst) - LOG_SQRT2PI - tconst * betaJoint * betaJoint / 2;
		logProbAltPoints		= logJointPDF + logScaleFactor;

	}// end of if

	return logProbNullPoints + logProbAltPoints;
}

/// <summary>
/// Compute Hvalues
/// Havlue is approximated Mvalue
/// </summary>
void MetaSnp::computeHvalues() 
{ // Hvalue is approximated Mvalue
	if (!isFixedEffectsComputed_) 
	{
		computeFixedEffects();
	}

	if (!isHvaluesComputed_) 
	{
		for (int i = 0; i < nStudy_; i++) 
		{
			MetaSnp metaSnp("dummy_rsid");
			
			for (int j = 0; j < nStudy_; j++) 
			{
				if (i != j) 
				{
					metaSnp.addStudy(betas_.at(j), standardErrors_.at(j));
				}
			}

			double var1 = pow(standardErrors_.at(i), 2.0) + pow(metaSnp.getStandardErrorFixedEffects(), 2.0);
			double var0 = pow(standardErrors_.at(i), 2.0);
			double f1 = exp(-pow(betas_.at(i) - metaSnp.getBetaFixedEffects(), 2.0) / (2.0*var1));
			double f0 = exp(-pow(betas_.at(i), 2.0) / (2.0*var0));
			double hvalue = f1 / (f0 + f1);
			hvalues_.push_back(hvalue);
		}

		isHvaluesComputed_ = true;
	}// end of if

	return;
}

/// <summary>
/// Print P and H values
/// </summary>
void MetaSnp::printPvaluesAndHvalues() 
{
	computeFixedEffects(1.0);
	computeRandomEffects();
	computeHanEskin(1.0, 1.0);

	if (pvalueHanEskinTabulated_ < pvalueFixedEffects_ && pvalueHanEskinTabulated_ < pvalueRandomEffects_) 
	{	
		printf("%s ", rsid_.c_str());

		for (int i = 0; i < nStudy_; i++) 
		{
			printf("%.5E ", getPvalue(i));
		}
		
		for (int i = 0; i < nStudy_; i++) 
		{
			printf("%.5f ", getHvalue(i));
		}

		printf("\n");
	}// end of if

	return;
}

/// <summary>
/// Compute Binary Effects
/// </summary>
void MetaSnp::computeBinaryEffectsStatistic() 
{
	if (!isHvaluesComputed_) 
	{
		computeHvalues();
	}
	
	double* zs = (double*)malloc(sizeof(double)*nStudy_);
	double* zWeights = (double*)malloc(sizeof(double)*nStudy_);

	for (int i = 0; i < nStudy_; i++) 
	{
		zs[i]		= betas_.at(i)  / standardErrors_.at(i);
		zWeights[i] = 1.0			/ standardErrors_.at(i);
	}

	double sumHvalueZweightZ = 0.0F;
	double sumHvalueSquareZweightSquare = 0.0F;

	for (int i = 0; i < nStudy_; i++) 
	{
		sumHvalueZweightZ += hvalues_.at(i) * zWeights[i] * zs[i];
		sumHvalueSquareZweightSquare += pow(hvalues_.at(i) * zWeights[i], 2.0F);
	}

	statisticBinaryEffects_ = sumHvalueZweightZ / sqrt(sumHvalueSquareZweightSquare);
	isBinaryEffectsStatisticComputed_ = true;

	free(zs);
	free(zWeights);

	return;
}

void MetaSnp::computeBinaryEffectsPvalue(long numSampling, unsigned int seed) 
{
	if (!isBinaryEffectsStatisticComputed_) 
	{
		computeBinaryEffectsStatistic();
	}
	
	srand(seed);
	double* zs = (double*)malloc(sizeof(double)*nStudy_);
	double* ws = (double*)malloc(sizeof(double)*nStudy_);

	for (int i = 0; i < nStudy_; i++) 
	{
		zs[i] = betas_.at(i) / standardErrors_.at(i);
		ws[i] = 1.0			 / standardErrors_.at(i);
	}

	double z = abs(statisticBinaryEffects_);
	// Code in this method currently follows C convention
	//		(i.e. lots of omissions in variable names! sorry.)
	int i, j, k, m, cnt, tmp;
	
	double pv, rnd, ratio, tau, lambda, ll, constz, ci95;
	double p_sample, p_situation, p_original, p_left, p_right;
	double mean, tmpsum, tmpsumsq, sumwssq;
	
	int n = nStudy_;
	
	double* xs			= (double*)malloc(sizeof(double)*(n + 1));
	double* stds		= (double*)malloc(sizeof(double)*(n + 1));
	double* pdf_num_alt = (double*)malloc(sizeof(double)*(n + 1));
	double* cdf_num_alt = (double*)malloc(sizeof(double)*(n + 1));
	double* samplezs	= (double*)malloc(sizeof(double)*(n + 1));
	double* samplebeta	= (double*)malloc(sizeof(double)*(n + 1));

	// Uniform probability of each scenario
	for (m = 1; m <= n; m++) 
	{
		pdf_num_alt[m] = 1.0 / n;
	}

	cdf_num_alt[1] = pdf_num_alt[1];

	for (m = 2; m <= n; m++) 
	{
		cdf_num_alt[m] = cdf_num_alt[m - 1] + pdf_num_alt[m];
	}

	// Get the factors of each scenario
	int* a = (int*)malloc(sizeof(int)*n);
	int x, y, t, done, rep;

	for (i = 0; i < n; i++) 
	{
		a[i] = i; 
	}

	for (m = 1; m <= n; m++) 
	{
		tmpsum = 0.0;
		tmpsumsq = 0.0;
		cnt = 0;

		for (rep = 0; rep < 1000; rep++) 
		{
			// permute "a" list
			for (i = 0; i < n; i++) 
			{
				j = i + (rand() % (n - i));
				std::swap(a[i], a[j]);
			}
			
			sumwssq = 0.0;
			
			for (i = 0; i < m; i++) 
			{
				sumwssq += ws[a[i]] * ws[a[i]];
			}
			
			for (i = 0; i < m; i++) 
			{
				mean = z * ws[a[i]] / sqrt(sumwssq);
				tmpsum += mean;
				tmpsumsq += mean * mean;
				cnt++;
			}

		}// end of for(rep)
		
		xs[m] = tmpsum / cnt;
		stds[m] = sqrt(tmpsumsq / cnt - xs[m] * xs[m] + 1);
	
	}
	
	// sample.
	pv = 0.0;
	cnt = 0;
	
	for (long sam = 0; sam < numSampling; sam++) 
	{
		// First, randomly choose situation(# of alt)
		rnd = rand() / (double)RAND_MAX;
		for (i = 1; i <= n; i++) 
		{
			if (rnd <= cdf_num_alt[i]) 
			{
				m = i; // "m" is the situation number.
				break;
			}
		}

		//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator(seed);
		std::normal_distribution<double> distribution(0.0, 1.0);

		// Given situation (m), sample n points
		for (i = 0; i < n; i++) 
		{
			samplezs[i] = distribution(generator);
		}

		for (i = 0; i < n; i++) 
		{
			if ((double)rand() / (double)RAND_MAX <= (double)m / n) 
			{
				samplezs[i] = samplezs[i] * stds[m] + xs[m];
			}
			samplebeta[i] = samplezs[i] / ws[i];
		}

		// run
		MetaSnp metaSnp("dummy_rsid");
		
		for (i = 0; i < n; i++) 
		{
			metaSnp.addStudy(samplebeta[i], 1.0 / ws[i]);
		}
		
		double sampleStatistic = metaSnp.getStatisticBinaryEffects();

		if (sampleStatistic >= z) 
		{ // significant.
			p_sample = 0.0; // pdf of sample distribution
			for (m = 1; m <= n; m++) 
			{
				p_situation = pdf_num_alt[m];
				
				for (k = 0; k < n; k++) 
				{
					p_left		= ((double)(n - m) / n)*(1. / sqrt(2 * M_PI))*exp(-0.5*samplezs[k] * samplezs[k]);
					p_right		= ((double)m / n)*(1. /sqrt(2 * M_PI*stds[m] * stds[m]))
								*exp(-0.5*(samplezs[k] - xs[m])*(samplezs[k] - xs[m]) / (stds[m] * stds[m]));
					p_situation *= p_left + p_right;
				}

				p_sample += p_situation;
			}
			
			p_original = 1.0;	// pdf of original distribution
			
			for (k = 0; k < n; k++)
			{ 
				p_original *= (1. / sqrt(2 * M_PI))*exp(-0.5*samplezs[k] * samplezs[k]);
			}
			
			ratio = p_original / p_sample;
			pv += ratio;
			cnt++;
		}// end of if
	}// end of for(sam)

	pv /= (double)numSampling;
	pv *= 2.;
	pvalueBinaryEffects_ = pv;
	
	if (pvalueBinaryEffects_ > 1.0) 
	{
		pvalueBinaryEffects_ = 1.0;
	}
	
	isBinaryEffectsPvalueComputed_ = true;

	free(xs);
	free(stds);
	free(pdf_num_alt);
	free(cdf_num_alt);
	free(samplezs);
	free(samplebeta);

	return;
}

std::string MetaSnp::configToString(bool* H1, int H1_n) 
{
	std::string str = "";

	for (int i = 0; i < H1_n; i++) 
	{
		str += H1[i] ? '1' : '0';
	}
	
	return str;
}

void MetaSnp::computeHanEskin(double lambdaMeanEffect, double lambdaHeterogeneity) 
{
	if (!isHeterogeneityComputed_) 
	{
		computeHeterogeneity();
	}

	if (!isPvalueTableRead_) 
	{
		printf("ERROR: Cannot compute HanEskin method without initializing p-value Table");
		exit(-1);
	}

	for (int i = 0; i < nStudy_; i++) 
	{
		betas[i]		= betas_.at(i);
		variances[i]	= pow(standardErrors_.at(i), 2.0);
		weights[i]		= 1.0 / variances[i];
	}

	double sumBetaProdWeight	= 0.0;
	double sumWeight			= 0.0;
	double sumWeightSquare		= 0.0;

	for (int i = 0; i < nStudy_; i++) 
	{
		sumBetaProdWeight	+= betas[i] * weights[i];
		sumWeight			+= weights[i];
		sumWeightSquare		+= weights[i] * weights[i];
	}
	
	// Iteratively find maximum likelihood parameters
	double previousMLBeta;
	double previousMLTauSquare;
	double previousLogLikelihood;
	double fixedEffectsMLBeta = sumBetaProdWeight / sumWeight;
	double nextMLBeta = fixedEffectsMLBeta;  // starting point
	double nextMLTauSquare = statisticTauSquare_; // starting point
	double nextLogLikelihood = NEGINF;//Double.NEGATIVE_INFINITY;
	double changeRatioMLBeta = POSINF;
	double changeRatioMLTauSquare = POSINF;
	double changeLogLikelihood = POSINF;
	double sumNumerator;
	double sumDenominator;

	while (changeRatioMLBeta > ML_ESTIMATE_CHANGE_RATIO_THRESHOLD || changeRatioMLTauSquare > ML_ESTIMATE_CHANGE_RATIO_THRESHOLD) 
	{
		previousMLBeta			= nextMLBeta;
		previousMLTauSquare		= nextMLTauSquare;
		previousLogLikelihood	= nextLogLikelihood;
		sumNumerator			= 0.0;
		sumDenominator			= 0.0;

		for (int i = 0; i < nStudy_; i++) 
		{
			sumNumerator	+= betas[i] / (variances[i] + previousMLTauSquare);
			sumDenominator	+=		1.0 / (variances[i] + previousMLTauSquare);
		}
		
		nextMLBeta		= sumNumerator / sumDenominator;
		sumNumerator	= 0.0;
		sumDenominator	= 0.0;

		for (int i = 0; i < nStudy_; i++) 
		{
			sumNumerator	+= (pow(betas[i] - nextMLBeta, 2.0) - variances[i]) / pow(variances[i] + previousMLTauSquare, 2.0);
			sumDenominator	+= 1.0 / pow(variances[i] + previousMLTauSquare, 2.0);
		}

		nextMLTauSquare = sumNumerator / sumDenominator;
		if (nextMLTauSquare < 0.0)
		{
			nextMLTauSquare = 0.0;
		}

		double ll = 0.0;
		
		for (int i = 0; i < nStudy_; i++) 
		{
			ll += 0.5 * log(2 * M_PI * (variances[i] + nextMLTauSquare)) - pow(betas[i] - nextMLBeta, 2.0) / (2 * (variances[i] + nextMLTauSquare));
		}
		
		nextLogLikelihood = ll;
		changeLogLikelihood = ll;
		changeLogLikelihood = nextLogLikelihood - previousLogLikelihood;
		changeRatioMLBeta = abs((nextMLBeta - previousMLBeta) / previousMLBeta);
		changeRatioMLTauSquare = abs((nextMLTauSquare - previousMLTauSquare) / previousMLTauSquare);

		if (changeLogLikelihood < 0.0) 
		{ // If somehow likelihood decreases,
			nextMLBeta		= previousMLBeta;  // step back and finish.
			nextMLTauSquare = previousMLTauSquare;
			break;
		}
	} // end of while

	double MLBeta		= nextMLBeta;
	double MLTauSquare	= nextMLTauSquare;
	// Compute statistics based on ML parameters
	double sumFormula1 = 0.0;
	double sumFormula2 = 0.0;
	double sumFormula3 = 0.0;
	double sumFormula4 = 0.0;

	for (int i = 0; i < nStudy_; i++) 
	{
		sumFormula1 += log(variances[i] / (variances[i] + MLTauSquare));
		sumFormula2 += betas[i] * betas[i] / variances[i];
		sumFormula3 += pow(betas[i] - fixedEffectsMLBeta, 2.0) / variances[i];
		sumFormula4 += pow(betas[i] - MLBeta, 2.0) / (variances[i] + MLTauSquare);
	}

	statisticHanEskinMeanEffectPart_	= sumFormula2 - sumFormula3;
	statisticHanEskinHeterogeneityPart_ = sumFormula1 + sumFormula3 - sumFormula4;
	
	if (statisticHanEskinHeterogeneityPart_ < 0.0) 
	{
		statisticHanEskinHeterogeneityPart_ = 0.0;
	}

	// Genomic-control
	statisticHanEskinMeanEffectPart_	/= lambdaMeanEffect;
	statisticHanEskinHeterogeneityPart_ /= lambdaHeterogeneity;
	statisticHanEskin_ = statisticHanEskinMeanEffectPart_ +	statisticHanEskinHeterogeneityPart_;

	// Compute asymptotic p-value
	pvalueHanEskinAsymptotic_ =
		0.5 * chiSquareComplemented(1.0, statisticHanEskin_) +//(1-boost::math::gamma_p<double, double>(1.0/2.0, statisticHanEskin_/2.0))+
		0.5 * chiSquareComplemented(2.0, statisticHanEskin_);//(1-boost::math::gamma_p<double, double>(2.0/2.0, statisticHanEskin_/2.0));
	
	// Use table to calculate accurate p-value
	if (nStudy_ <= TABLE_MAX_NSTUDY) 
	{
		int nearestIndexBottom = floor(statisticHanEskin_ * 10.0);
		if (nearestIndexBottom < 0)
		{
			nearestIndexBottom = 0;
		}

		int nearestIndexTop = ceil(statisticHanEskin_ * 10.0);
		if (nearestIndexTop < 0)
		{
			nearestIndexTop = 0;
		}

		if (nearestIndexTop < TABLE_NCOLUMN) 
		{
			int rowNumber = nStudy_ - 2;
			double tablePvalueAtIndexBottom = 0.0;
			
			try 
			{
				tablePvalueAtIndexBottom = pvalueTable_[rowNumber][nearestIndexBottom];
			}
			catch (std::exception e) 
			{
				printf("%f %f %d %d\n", 
					statisticHanEskinMeanEffectPart_,
					statisticHanEskinHeterogeneityPart_,
					nearestIndexBottom,
					nearestIndexTop);
				exit(-1);
			}

			double asymptoticPvalueAtIndexBottom =
				0.5 * chiSquareComplemented(1.0, nearestIndexBottom / 10.0) +//(1-boost::math::gamma_p<double, double>(1.0/2.0, (nearestIndexBottom / 10.0)/2.0)) +
				0.5 * chiSquareComplemented(2.0, nearestIndexBottom / 10.0);//(1-boost::math::gamma_p<double, double>(2.0/2.0, (nearestIndexBottom / 10.0)/2.0));
			
			double ratioAtIndexBottom = tablePvalueAtIndexBottom /
										asymptoticPvalueAtIndexBottom;
			
			double tablePvalueAtIndexTop = pvalueTable_[rowNumber][nearestIndexTop];
			
			double asymptoticPvalueAtIndexTop =
				0.5 * chiSquareComplemented(1.0, nearestIndexTop / 10.0) +//(1-boost::math::gamma_p<double, double>(1.0/2.0, (nearestIndexTop / 10.0)/2.0)) +
				0.5 * chiSquareComplemented(2.0, nearestIndexTop / 10.0);//(1-boost::math::gamma_p<double, double>(2.0/2.0, (nearestIndexTop / 10.0)/2.0));
			
			double ratioAtIndexTop = tablePvalueAtIndexTop /
										asymptoticPvalueAtIndexTop;
			
			double ratioInterpolated =
				ratioAtIndexBottom + (ratioAtIndexTop - ratioAtIndexBottom) *
				(statisticHanEskin_ - nearestIndexBottom / 10.0) / 0.1;
			
			pvalueHanEskinTabulated_ = ratioInterpolated * pvalueHanEskinAsymptotic_;
		}
		else 
		{
			int    rowNumber = nStudy_ - 2;
			double tablePvalueAtTheEnd = pvalueTable_[rowNumber][TABLE_NCOLUMN - 1];
		
			double asymptoticPvalueAtTheEnd =
				0.5 * chiSquareComplemented(1.0, TABLE_MAX_THRESHOLD) +//(1-boost::math::gamma_p<double, double>(1.0/2.0, TABLE_MAX_THRESHOLD/2.0)) +
				0.5 * chiSquareComplemented(2.0, TABLE_MAX_THRESHOLD);//(1-boost::math::gamma_p<double, double>(2.0/2.0, TABLE_MAX_THRESHOLD/2.0));
		
			double ratioAtTheEnd = tablePvalueAtTheEnd /
									asymptoticPvalueAtTheEnd;
		
			pvalueHanEskinTabulated_ = ratioAtTheEnd * pvalueHanEskinAsymptotic_;
		}
	}
	else 
	{
		pvalueHanEskinTabulated_ = pvalueHanEskinAsymptotic_;
	}
	
	isHanEskinComputed_ = true;

	return;
}

void MetaSnp::computeHanEskin() 
{
	computeHanEskin(1.0, 1.0);
	return;
}

double MetaSnp::getHvalue(int i) 
{
	if (!isHvaluesComputed_)
	{
		computeHvalues();
	}

	return hvalues_.at(i);
}

double MetaSnp::getMvalue(int i) 
{
	if (!isMvaluesComputed_) 
	{
		// computeMvalues();
		printf("First compute m-value before calling");
		exit(-1);
	}

	return mvalues_.at(i);
}

double MetaSnp::getPvalue(int i) 
{
	return chiSquareComplemented(1.0, pow(betas_.at(i) / standardErrors_.at(i), 2.0));
}

double MetaSnp::getPvalueFixedEffects() 
{
	if (!isFixedEffectsComputed_) 
	{
		computeFixedEffects();
	}
	return pvalueFixedEffects_;
}

double MetaSnp::getBetaFixedEffects() 
{
	if (!isFixedEffectsComputed_) 
	{
		computeFixedEffects();
	}
	
	return betaFixedEffects_;
}

double MetaSnp::getStandardErrorFixedEffects() 
{
	if (!isFixedEffectsComputed_) 
	{
		computeFixedEffects();
	}
	
	return standardErrorFixedEffects_;
}

double MetaSnp::getPvalueRandomEffects() 
{
	if (!isRandomEffectsComputed_) 
	{
		computeRandomEffects();
	}
	
	return pvalueRandomEffects_;
}

double MetaSnp::getBetaRandomEffects() 
{
	if (!isRandomEffectsComputed_) 
	{
		computeRandomEffects();
	}
	
	return betaRandomEffects_;
}

double MetaSnp::getStandardEfforRandomEffects() 
{
	if (!isRandomEffectsComputed_) 
	{
		computeRandomEffects();
	}
	
	return standardErrorRandomEffects_;
}

double MetaSnp::getPvalueHanEskin() 
{
	if (!isHanEskinComputed_) 
	{
		computeHanEskin();
	}
	
	return pvalueHanEskinTabulated_;
}

double MetaSnp::getStatisticHanEskinMeanEffectPart() 
{
	if (!isHanEskinComputed_) 
	{
		computeHanEskin();
	}
	
	return statisticHanEskinMeanEffectPart_;
}

double MetaSnp::getStatisticHanEskinHeterogeneityPart() 
{
	if (!isHanEskinComputed_) 
	{
		computeHanEskin();
	}
	
	return statisticHanEskinHeterogeneityPart_;
}

double MetaSnp::getStatisticBinaryEffects() 
{
	if (!isHanEskinComputed_) 
	{
		computeHanEskin();
	}
	
	return statisticHanEskinHeterogeneityPart_;
}

double MetaSnp::getPvalueBinaryEffects() 
{
	if (!isBinaryEffectsPvalueComputed_) 
	{
		printf("First compute binary effects p-value before calling");
		exit(-1);
	}
	
	return pvalueBinaryEffects_;
}

void MetaSnp::printHeadings(FILE* outFile) 
{
	std::string str = "RSID\t#STUDY\tPVALUE_FE\tBETA_FE\tSTD_FE\tPVALUE_RE\tBETA_RE\tSTD_RE\t";
	str += "PVALUE_RE2\tSTAT1_RE2\tSTAT2_RE2\tPVALUE_BE\tI_SQUARE\tQ\tPVALUE_Q\tTAU_SQUARE\t";
	str += "PVALUES_OF_STUDIES(Tab_delimitered)\tMVALUES_OF_STUDIES(Tab_delimitered)\n";
	fprintf(outFile,"%s",str.c_str());

	return;
}

void MetaSnp::printResults(FILE* f) 
{
	fprintf(f,"%s\t", rsid_.c_str());
	fprintf(f,"%d\t", nStudy_);
	
	if (isFixedEffectsComputed_) 
	{
		fprintf(f,"%.6G\t", pvalueFixedEffects_);
		fprintf(f,"%.6G\t", betaFixedEffects_);
		fprintf(f,"%.6G\t", standardErrorFixedEffects_);
		fprintf(f,"%.6G\t", pvalueRandomEffects_);
		fprintf(f,"%.6G\t", betaRandomEffects_);
		fprintf(f,"%.6G\t", standardErrorRandomEffects_);
		fprintf(f,"%.6G\t", pvalueHanEskinTabulated_);
		fprintf(f,"%.6G\t", statisticHanEskinMeanEffectPart_);
		fprintf(f,"%.6G\t", statisticHanEskinHeterogeneityPart_);
		
		if (isBinaryEffectsPvalueComputed_) 
		{
			fprintf(f,"%.6G\t", pvalueBinaryEffects_);
		}
		else 
		{
			fprintf(f,"NA\t");
		}
		
		fprintf(f,"%.6G\t", statisticISquare_);
		fprintf(f,"%.6G\t", statisticQ_);
		fprintf(f,"%.6G\t", pvalueQ_);
		fprintf(f,"%.6G\t", statisticTauSquare_);
	}
	else 
	{ // if not, it must be a problematic SNPs with nStudy < 2; just print NA
		for (int i = 0; i < 14; i++) 
		{
			fprintf(f,"NA\t");
		}
	}

	int j;
	j = 0;
	for (int i = 0; i < nStudyIncludingNa_; i++) 
	{
		if (isNa_.at(i)) 
		{
			fprintf(f,"NA\t");
		}
		else 
		{
			fprintf(f,"%G\t", getPvalue(j));
			++j;
		}
	}

	j = 0;
	for (int i = 0; i < nStudyIncludingNa_; i++) 
	{
		if (isNa_.at(i))
		{
			fprintf(f,"NA\t");
		}
		else 
		{
			if (isMvaluesComputed_) 
			{
				fprintf(f,"%.3f\t", getMvalue(j));
			}
			else 
			{
				fprintf(f,"NA\t");
			}
			++j;
		}// end of if-else
	}

	fprintf(f, "\n");

	return;
}

void MetaSnp::readPvalueTableFile(std::string pvalueTableFile) 
{
	pvalueTable_ = (double**)malloc(sizeof(double*)*TABLE_NROW);
	
	for (int i = 0; i < TABLE_NROW; i++) 
	{
		pvalueTable_[i] = (double*)malloc(sizeof(double)*TABLE_NCOLUMN);
	}
	
	std::ifstream infile;
	try 
	{
		infile = std::ifstream(pvalueTableFile);
	}
	catch (std::exception e) 
	{
		std::cout << "ERROR: P-value Table file Error" << std::endl;
		exit(-1);
	}

	try 
	{
		std::string readLine;
		std::getline(infile, readLine); // ignore top row
		std::vector<std::string> tokens;
		tokens.resize(TABLE_NCOLUMN + 1);
	
		for (int i = 0; i < TABLE_NROW; i++) 
		{
			if (!std::getline(infile, readLine)) 
			{
				std::cout << "ERROR: Reading error from P-value Table file" << std::endl;
				exit(-1);
			}
			
			split(tokens, readLine, " ");
			
			if (tokens.size() < TABLE_NCOLUMN+1) 
			{
				std::cout << "ERROR: P-value Table File has too few columns" << std::endl;
				exit(-1);
			}

			for (int j = 0; j < TABLE_NCOLUMN; j++) 
			{
				try 
				{
					pvalueTable_[i][j] = std::stod(tokens.at(j+1));
				}
				catch (std::exception e) 
				{
					std::cout << "Incorrect float value in Pvalue Table file."<<std::endl;
					exit(-1);
				}
			} // end of for(j)
			tokens.clear();
		} // end of for(i)
	}
	catch (std::exception e) 
	{
		std::cout << "ERROR: error encountered while reading Pvalue Table file" << std::endl;
		exit(-1);
	}

	isPvalueTableRead_ = true;

	return;
}

// END OF FILE