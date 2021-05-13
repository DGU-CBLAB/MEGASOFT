#include"Scheduler.h"

VOID Scheduler::prepare() // prev MetaAnalysis Prep
{
	std::cout << "----- Preparing Analysis\n";

	prepareReadPvalueTableFile();
	
	srand(mp_option->getSeed());
	
	prepareInitVal();
	prepareOpenFile();
	prepareReadInputFile();

	ResultPrinter::printResultHeadings(mp_outFile);

	return;
}
VOID Scheduler::process()
{
	doMetaAnalysis();
	computeLambda();
}
VOID Scheduler::postProcess()
{
	postCloseFile();
	postReorderOutput();
}
// Prepare
VOID Scheduler::prepareReadPvalueTableFile()
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


	try {
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
VOID Scheduler::prepareInitVal()
{
	mp_option->setNumSnps(0);
	mp_option->setMaxNumStudy(0);

	mvd_heterogeneityParts.clear();
	mvd_meanEffectParts.clear();
}
VOID Scheduler::prepareOpenFile()
{
	try
	{
		//mp_inFile = fopen(mp_option->getInputFile().c_str(), FILE_READ);
		mp_outFile = fopen(mp_option->getOutputFile().c_str(), FILE_WRITE);
	}
	catch (std::exception e)
	{
		Options::printErrorAndQuit(e.what());
	}
}
VOID Scheduler::prepareReadInputFile()
{
	try
	{
		std::ifstream inStream(mp_option->getInputFile());

		UINT32 doneCnt = 0;
		std::string line;

		while (std::getline(inStream, line))
		{
			mvs_inputLines.push_back(line);
		}

		inStream.close();
	}
	catch (std::exception e)
	{
		Options::printErrorAndQuit(e.what());
	}
}

// Process
VOID Scheduler::doMetaAnalysis()
{
	std::cout << "----- Performing meta-analysis\n";

}
VOID Scheduler::computeLambda()
{
	std::cout << "---- Performing lambda compute\n";
	DOUBLE median;
	DOUBLE expectedMedian;
	if (mvd_meanEffectParts.size() > 0)
	{
		std::sort(mvd_meanEffectParts.begin(), mvd_meanEffectParts.end());
		median = mvd_meanEffectParts.at((int)(mvd_meanEffectParts.size() / 2.0));
		expectedMedian = pow(boost::math::find_location<boost::math::normal>(boost::math::complement(0, 0.25, 1.0)), 2.0);
		mp_option->setOutputLambdaMeanEffect(median / expectedMedian);
	}
	if (mvd_heterogeneityParts.size() > 0)
	{
		std::sort(mvd_heterogeneityParts.begin(), mvd_heterogeneityParts.end());
		median = mvd_heterogeneityParts.at((int)(mvd_heterogeneityParts.size() / 2.0));
		if (mp_option->getMaxNumStudy() > 50)
		{
			expectedMedian = pow(boost::math::find_location<boost::math::normal>(boost::math::complement(0, 0.25, 1.0)), 2.0);
		}
		else
		{
			expectedMedian = expectedMedianHanEskinHeterogeneityPart[mp_option->getMaxNumStudy() - 2];
		}
		mp_option->setOutputLambdaHeterogeneity(median / expectedMedian);
	}
	return;
}

// Post
VOID Scheduler::postCloseFile()
{
	try
	{
		//fclose(mp_inFile);
		fclose(mp_outFile);
	}
	catch (std::exception e)
	{
		Options::printErrorAndQuit(e.what());
	}
}
VOID Scheduler::postReorderOutput()
{
	try {
		std::string tmp_str;
		FILE* file = fopen(mp_option->getOutputFile().c_str(), "r");
		std::ifstream Instream(mp_option->getOutputFile().c_str());

		std::getline(Instream, tmp_str); // ignore first line

		std::vector<ThreadResult> total;

		while (std::getline(Instream, tmp_str))
		{
			ThreadResult mt(stoi(tmp_str.substr(0, tmp_str.find('\t'))), tmp_str);
			total.push_back(mt);
		}

		sort(total.begin(), total.end(), compareKeyValues);
		Instream.close();
		FILE* outfile = fopen(outputFile_.c_str(), "w");
		MetaSnp::printHeadings(outfile);
		for (int i = 0; i < total.size(); i++)
		{
			fprintf(outfile, "%s\n", total.at(i).values_str.c_str());
		}

		fclose(outfile);
	}
	catch (std::exception e)
	{
		printf("ERROR: Posterior.txt file has been altered!\n[%s]", e.what());
		std::exit(-1);
	}
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