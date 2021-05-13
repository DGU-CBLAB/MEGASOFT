#pragma once
#include"Options.h"
#include"CPU_Analysis.h"
#include<boost/math/distributions/find_location.hpp>
#include<boost/math/distributions/normal.hpp>
#include<algorithm>
#include<fstream>
#include<random>

static const double expectedMedianHanEskinHeterogeneityPart[] = // from nStudy 2 to 50
{	0.2195907137,0.2471516439,0.2642270318,0.2780769264,0.2886280267,
	0.2977812664,0.3020051148,0.3091428179,0.3158605559,0.3221788173,
	0.3259133140,0.3295976587,0.3335375196,0.3358395088,0.3368309971,
	0.3421941686,0.3448030927,0.3463590948,0.3477384754,0.3487900288,
	0.3494938171,0.3542087791,0.3573286353,0.3589703411,0.3586951356,
	0.3596101209,0.3605611682,0.3624799993,0.3648322669,0.3659817739,
	0.3671267389,0.3693952373,0.3693395144,0.3696863113,0.3706067524,
	0.3718103285,0.3749536619,0.3758886239,0.3753612342,0.3781458299,
	0.3798346038,0.3763434983,0.3796968747,0.3784334922,0.3794411347,
	0.3808582942,0.3813485882,0.3843230993,0.3824863479 
};

class Scheduler
{
private:
	Options* mp_option;

	//FILE* mp_inFile;
	FILE* mp_outFile;
	std::vector<DOUBLE> mvd_meanEffectParts;
	std::vector<DOUBLE> mvd_heterogeneityParts;
	std::vector<std::string> mvs_inputLines;

	DOUBLE ma_pvalueTable[HANESKIN_TABLE_ROW][HANESKIN_TABLE_COL];

private:
	VOID split(std::vector<std::string>& tokens, const std::string& str, const std::string& delim);
// Prepare
	VOID prepareReadPvalueTableFile();
	VOID prepareInitVal();
	VOID prepareOpenFile();
	VOID prepareReadInputFile();
// Process
	VOID doMetaAnalysis();
	VOID computeLambda();

// Post
	VOID postCloseFile();
	VOID postReorderOutput();

public:
	Scheduler(Options* option)
		:mp_option(option) {}
	VOID prepare();
	VOID process();
	VOID postProcess();

};