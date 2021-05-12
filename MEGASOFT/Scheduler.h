#pragma once
#include"Options.h"
class Scheduler
{
private:
	Options* mp_option;
	DOUBLE ma_pvalueTable[HANESKIN_TABLE_ROW][HANESKIN_TABLE_COL];

	VOID readPvalueTableFile();
	VOID doMetaAnalysis();
	VOID computeLambda();

	VOID split(std::vector<std::string>& tokens, const std::string& str, const std::string& delim);
public:
	Scheduler(Options* option)
		:mp_option(option) {}
	VOID prepare();
	VOID run();
	VOID cleanUp();
};