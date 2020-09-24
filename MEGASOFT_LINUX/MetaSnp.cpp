#include"MetaSnp.h"
#include<boost/math/special_functions/beta.hpp>
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

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


void MetaSnp::readPvalueTableFile(std::string pvalueTableFile)
{
    try{
        std::ifstream infile(pvalueTableFile);

        std::string readLine;
        std::getline(infile, readLine); // ignore top line
        
        std::vector<std::string> tokens;
        tokens.resize(TABLE_NCOLUMN+1);

        for(int r =0;r<TABLE_NROW; r++)
        {
            std::getline(infile, readLine);

            split(tokens, readLine, " ");

            for(int c=0;c<TABLE_NCOLUMN;c++)
            {
                pvalueTable_[r][c] = std::stod(tokens.at(c+1));

            }
            tokens.clear();
        }

    }
    catch(ifstream::failure e)
    {
        printf("Error: P-value Table file Error [%s]\n", pvalueTableFile);
        exit(ERROR_IO);
    }
    catch(std::exception e){
        printf("Error: Un recognizable Pvalue Table File [%s]\n", pvalueTableFile);
        exit(ERROR_PVALUE_FORMAT_ERROR);
    }

    isPvalueTableRead_ = true;

    return;
}