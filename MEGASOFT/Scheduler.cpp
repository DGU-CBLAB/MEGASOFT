#include"Scheduler.h"
#include<fstream>
VOID Scheduler::prepare()
{
	readPvalueTableFile();
}
VOID Scheduler::run()
{
	doMetaAnalysis();
	computeLambda();
}
VOID Scheduler::cleanUp()
{

}

VOID Scheduler::readPvalueTableFile() 
{
	std::ifstream infile;
	std::string readLine;
	std::vector<std::string> tokens;

	try 
	{
		infile = std::ifstream(mp_option->getPvalueTableFile());
	}
	catch (std::exception e)
	{
		std::cout << "ERROR: P-value Table file Error" << std::endl;
		exit(ERROR_FILE_NOT_FOUND);
	}


	try{
		std::getline(infile, readLine); // ignore top row
		tokens.resize(HANESKIN_TABLE_COL + 1);

		for (int row = 0; row < HANESKIN_TABLE_ROW; row++)
		{
			if (!std::getline(infile, readLine))
			{
				std::cout << "ERROR: Reading error from P-value Table file" << std::endl;
				exit(ERROR_FILE_READ_FAIL);
			}

			split(tokens, readLine, " ");

			if (tokens.size() < HANESKIN_TABLE_COL + 1)
			{
				std::cout << "ERROR: P-value Table File has too few columns" << std::endl;
				exit(ERROR_FILE_READ_FAIL);
			}

			for (int col = 0; col < HANESKIN_TABLE_COL; col++)
			{
				try
				{
					ma_pvalueTable[row][col] = stod(tokens.at(col + 1));
				}
				catch (std::exception e)
				{
					std::cout << "Incorrect float value in Pvalue Table file." << std::endl;
					exit(ERROR_FILE_READ_FAIL);
				}
			}
			tokens.clear();
		}
	}
	catch (std::exception e)
	{
		std::cout << "ERROR: error encountered while reading Pvalue Table file" << std::endl;
		exit(ERROR_FILE_READ_FAIL);
	}
	return;
}
VOID Scheduler::doMetaAnalysis()
{
	std::cout << "----- Performing meta-analysis\n";
}
VOID Scheduler::computeLambda()
{
	std::cout << "---- Performing lambda compute\n";
}

VOID Scheduler::split(std::vector<std::string>& tokens, const std::string& str, const std::string& delim)
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