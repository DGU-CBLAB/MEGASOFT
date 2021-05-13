#include"CPU_Analysis.h"

CPU_Analysis::CPU_Analysis(std::string line)
{
	split(mvs_tokens, line, DELIM_SPACE);

}
VOID CPU_Analysis::process()
{
	if (mvs_tokens.size() > 1
		|| mvs_tokens.at(0).at(0) != CHAR_SHARP)
	{
		return;
	}


	return;
}
VOID CPU_Analysis::postProcess(FILE* outFile)
{
	ResultPrinter::printResults(outFile, &m_analysisData);
}
VOID CPU_Analysis::split(std::vector<std::string>& tokens, const std::string& str, const std::string& delim)
{
	size_t prev = 0, pos = 0;
	// Empty token results
	tokens.clear();

	do
	{
		pos = str.find(delim, prev);

		if (str.find(DELIM_TAB, prev) < pos)
		{
			pos = str.find(DELIM_TAB, prev);
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