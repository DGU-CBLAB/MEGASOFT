#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "Gamma.cu"
#define M_PI acos(-1.0)

static const int TABLE_NROW = 49;
static const int TABLE_MAX_NSTUDY = 50;
static const int TABLE_NCOLUMN = 331;
static double TABLE_MAX_THRESHOLD;
static double** pvalueTable_;
static bool isPvalueTableRead_; // TABLE_NROW x TABLE_NCOLUMN

typedef struct Options {
	int nStudy_ = 0;
	int nStudyIncludingNa_ = 0;
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

	static double ML_ESTIMATE_CHANGE_RATIO_THRESHOLD;
	static double LOG_SQRT2PI;
}Options;
typedef struct Data {
	double betas_;
	double standardErrors_;
	bool   isNa_;
	double hvalues_;
	double mvalues_;
	double betas, variances, weights;
}Data;

__device__ __host__ void gpu_logBeta(double* res, double m, double n);
__device__ __host__ void gpu_chiSquareComplemented(double* res, double v, double x);