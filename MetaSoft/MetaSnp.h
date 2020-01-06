#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include<iostream>
#include<stdio.h>
#include<sstream>
#include<fstream>
#include<stdio.h>
#include<string>
#include<chrono>
#include<random>
#include<fstream>
#include<vector>
#include<time.h>
#include<math.h>
#include<cmath>
#include<mutex>
#include<thread>
#include<map>
#include<boost/math/distributions/normal.hpp>
#include<boost/math/distributions/beta.hpp>
#include<boost/math/distributions/find_location.hpp>
#include<boost/math/distributions.hpp>
//#define M_PI acos(-1.0) // Accurate PI constant
#define NORMAL_EXECUTION 1
#define ABNORMAL_EXECUTION -1
// Statistical functions

class map_tuple {
public:
	int key;
	std::string val;
	map_tuple(int key, std::string val) {
		this->key = key;
		this->val = val;
	}
};
bool map_comp(const map_tuple& a, const map_tuple& b);
void split(std::vector<std::string>& tokens, const std::string& str, const std::string& delim);

class MetaSnp {
private:
	std::string  rsid_;
	int nStudy_				= 0;
	int nStudyIncludingNa_	= 0;
	double statisticFixedEffects_;
	double pvalueFixedEffects_;
	double betaFixedEffects_;
	double standardErrorFixedEffects_;
	double statisticRandomEffects_;
	double pvalueRandomEffects_;
	double betaRandomEffects_;
	double standardErrorRandomEffects_;
	double statisticHanEskin_;
	double statisticHanEskinMeanEffectPart_;
	double statisticHanEskinHeterogeneityPart_;
	double pvalueHanEskinAsymptotic_;
	double pvalueHanEskinTabulated_;
	double statisticBinaryEffects_;
	double pvalueBinaryEffects_;
	double statisticQ_;
	double pvalueQ_;
	double statisticISquare_;
	double statisticTauSquare_;
	bool isFixedEffectsComputed_ = false;
	bool isRandomEffectsComputed_ = false;
	bool isHeterogeneityComputed_ = false;
	bool isHvaluesComputed_ = false;
	bool isMvaluesComputed_ = false;
	bool isHanEskinComputed_ = false;
	bool isBinaryEffectsStatisticComputed_ = false;
	bool isBinaryEffectsPvalueComputed_ = false;
	std::vector<double>  betas_;
	std::vector<double>  standardErrors_;
	std::vector<bool>	 isNa_;
	std::vector<double>  hvalues_;
	std::vector<double>  mvalues_;
	static double ML_ESTIMATE_CHANGE_RATIO_THRESHOLD;
	static double LOG_SQRT2PI;
public:
	double logBeta(double m, double n);
	double chiSquareComplemented(double x, double v);
	MetaSnp(std::string rsid);
	void addStudy(double beta, double standardError);
	void addNaStudy();
	void computeFixedEffects(double lambdaMeanEffect);
	void computeFixedEffects();
	void computeHeterogeneity();
	void computeRandomEffects();
	void computeMvalues(double priorAlpha, double priorBeta, double priorSigma);
	void computeMvaluesMCMC(double priorAlpha, double priorBeta, double priorSigma,
		long sample, long burnin, double probRandom, double maxNumFlipArg, int seed);
private:
	double observationLogLikelihood(double* betas, int betas_size, double* ts, bool* H1, int numH1, double priorVar);
public:
	void computeHvalues();
	void printPvaluesAndHvalues();
	void computeBinaryEffectsStatistic();
	void computeBinaryEffectsPvalue(long numSampling, int seed);
	void computeHanEskin(double lambdaMeanEffect, double lambdaHeterogeneity);
protected:
	static std::string configToString(bool* H1, int H1_n);
public:
	void computeHanEskin();
	int getNStudy() { return nStudy_; }
	double getHvalue(int i);
	double getMvalue(int i);
	double getPvalue(int i);
	double getPvalueFixedEffects();
	double getBetaFixedEffects();
	double getStandardErrorFixedEffects();
	double getPvalueRandomEffects();
	double getBetaRandomEffects();
	double getStandardEfforRandomEffects();
	double getPvalueHanEskin();
	double getStatisticHanEskinMeanEffectPart();
	double getStatisticHanEskinHeterogeneityPart();
	double getStatisticBinaryEffects();
	double getPvalueBinaryEffects();
	static void printHeadings(FILE* outFile);
	void printResults(FILE* f); // used to be static method

public:
	// Class-variable part for
	// P-value Table of Han Eskin Statistic.
	// Currently, table dimension is fixed as follows.
	// nStudy (rows) from 2 to 50 (toal 49 rows)
	// statistic thresholds (columns) from 0.0 to 33.0 (total 331 columns)
	static const int TABLE_NROW = 49;
	static const int TABLE_MAX_NSTUDY = 50;
	static const int TABLE_NCOLUMN = 331;
	static double TABLE_MAX_THRESHOLD;
	static double** pvalueTable_;
	static bool isPvalueTableRead_; // TABLE_NROW x TABLE_NCOLUMN

public:
	static void readPvalueTableFile(std::string pvalueTableFile);
};