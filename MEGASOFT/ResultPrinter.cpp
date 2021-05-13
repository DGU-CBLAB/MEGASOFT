#pragma once
#include"ResultPrinter.h"

VOID ResultPrinter::printResultHeadings(FILE* outFile)
{
	std::string str = "RSID\t#STUDY\tPVALUE_FE\tBETA_FE\tSTD_FE\tPVALUE_RE\tBETA_RE\tSTD_RE\t";
	str += "PVALUE_RE2\tSTAT1_RE2\tSTAT2_RE2\tPVALUE_BE\tI_SQUARE\tQ\tPVALUE_Q\tTAU_SQUARE\t";
	str += "PVALUES_OF_STUDIES(Tab_delimitered)\tMVALUES_OF_STUDIES(Tab_delimitered)\n";
	fprintf(outFile, "%s", str.c_str());
}
VOID ResultPrinter::printResults(FILE* f, AnalysisData* p_analysisData)
{
	//fprintf(f, "%s\t", rsid_.c_str());
	//fprintf(f, "%d\t", nStudy_);

	//if (isFixedEffectsComputed_)
	//{
	//	fprintf(f, "%.6G\t", pvalueFixedEffects_);
	//	fprintf(f, "%.6G\t", betaFixedEffects_);
	//	fprintf(f, "%.6G\t", standardErrorFixedEffects_);
	//	fprintf(f, "%.6G\t", pvalueRandomEffects_);
	//	fprintf(f, "%.6G\t", betaRandomEffects_);
	//	fprintf(f, "%.6G\t", standardErrorRandomEffects_);
	//	fprintf(f, "%.6G\t", pvalueHanEskinTabulated_);
	//	fprintf(f, "%.6G\t", statisticHanEskinMeanEffectPart_);
	//	fprintf(f, "%.6G\t", statisticHanEskinHeterogeneityPart_);

	//	if (isBinaryEffectsPvalueComputed_)
	//	{
	//		fprintf(f, "%.6G\t", pvalueBinaryEffects_);
	//	}
	//	else
	//	{
	//		fprintf(f, "NA\t");
	//	}

	//	fprintf(f, "%.6G\t", statisticISquare_);
	//	fprintf(f, "%.6G\t", statisticQ_);
	//	fprintf(f, "%.6G\t", pvalueQ_);
	//	fprintf(f, "%.6G\t", statisticTauSquare_);
	//}
	//else
	//{ // if not, it must be a problematic SNPs with nStudy < 2; just print NA
	//	for (int i = 0; i < 14; i++)
	//	{
	//		fprintf(f, "NA\t");
	//	}
	//}

	//int j;
	//j = 0;
	//for (int i = 0; i < nStudyIncludingNa_; i++)
	//{
	//	if (isNa_.at(i))
	//	{
	//		fprintf(f, "NA\t");
	//	}
	//	else
	//	{
	//		fprintf(f, "%G\t", getPvalue(j));
	//		++j;
	//	}
	//}

	//j = 0;
	//for (int i = 0; i < nStudyIncludingNa_; i++)
	//{
	//	if (isNa_.at(i))
	//	{
	//		fprintf(f, "NA\t");
	//	}
	//	else
	//	{
	//		if (isMvaluesComputed_)
	//		{
	//			fprintf(f, "%.3f\t", getMvalue(j));
	//		}
	//		else
	//		{
	//			fprintf(f, "NA\t");
	//		}
	//		++j;
	//	}// end of if-else
	//}

	//fprintf(f, "\n");

	//return;

}
