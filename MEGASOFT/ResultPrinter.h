#pragma once
#include"Config.h"
#include"AnalysisData.h"
class ResultPrinter
{
private:

public:
	static VOID printResultHeadings(FILE* outFile);
	static VOID printResults(FILE* outFile, AnalysisData* p_option);
};