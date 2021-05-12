#include"Config.h"
#include"Log.h"
#include <boost/program_options.hpp>

class Options:public Log
{
private:
	std::string m_inputFile;
	std::string m_outputFile;
	std::string m_pvalueTableFile;
	std::string m_logFile;
	DOUBLE	m_inputLambdaMeanEffect;
	DOUBLE	m_inputLambdaHeterogeneity;
	BOOL	m_willComputeMvalue;
	DOUBLE	m_priorSigma;
	DOUBLE	m_priorAlpha;
	DOUBLE	m_priorBeta;
	DOUBLE	m_mvaluePvalueThreshold;
	std::string  m_mvalueMethod;
	LONG    m_mcmcSample;
	LONG    m_mcmcBurnin;
	DOUBLE	m_mcmcProbRandom;
	DOUBLE	m_mcmcMaxNumFlip;
	BOOL	m_willComputeBinaryEffects;
	LONG    m_binaryEffectsSample;
	LONG    m_binaryEffectsLargeSample;
	DOUBLE	m_binaryEffectsPvalueThreshold;
	UINT32	m_seed;
	BOOL	m_isVerbose;
	UINT32	m_CPUNumThread;
	BOOL	m_gpu;

	UINT32	m_numSnps;
	UINT32	m_maxNumStudy;
	DOUBLE	m_outputLambdaMeanEffect;
	DOUBLE	m_outputLambdaHeterogeneity;
public:
	Options()
	{
		initialize();
	}
	VOID initialize();
	VOID handleArguments(INT32 argc, char* argv[]);
	VOID createSummary(INT32 argc, char* argv[]);
	
	VOID printLog(DOUBLE time);
	VOID printArguments(VOID);
	VOID printErrorAndQuit(std::string msg);

	std::string getPvalueTableFile(VOID) { return m_pvalueTableFile; }
};