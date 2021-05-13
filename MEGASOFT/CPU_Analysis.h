#pragma once
#include"ResultPrinter.h"
#include<vector>
class CPU_Analysis
{
private:
	AnalysisData m_analysisData;
	std::vector<std::string> mvs_tokens;

	VOID split(std::vector<std::string>& tokens, const std::string& str, const std::string& delim);

public:
	CPU_Analysis(std::string line);

	VOID process(VOID);
	VOID postProcess(FILE* outFile);

};