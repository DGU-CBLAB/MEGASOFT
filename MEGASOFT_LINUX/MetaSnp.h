#include<iostream>
#include<stdio.h>
#include<sstream>
#include<fstream>
#include<stdio.h>
#include<string>
#include<chrono>
#include<random>
#include<fstream>
#include<vector>
#include<time.h>
#include<math.h>
#include<cmath>
#include<mutex>
#include<thread>
#include<pthread.h>
#include<map>

// Boost Library 1.72.1
#include<boost/math/distributions/normal.hpp>
#include<boost/math/distributions/beta.hpp>
#include<boost/math/distributions/find_location.hpp>
#include<boost/math/distributions.hpp>


#define THREAD (3)

#define ERROR_IO					(0)
#define ERROR_IO_FILE_OPEN_CLOSE	(1)
#define ERROR_TREHAD_CREATE 		(2)
#define ERROR_THREAD_JOIN			(3)
#define ERROR_META_ANALYSIS			(4)
#define ERROR_PVALUE_FORMAT_ERROR	(5)
#define ERROR_UNDEFINED				(6)

const char* ERROR_MESSAGE[]={
	"ERROR_IO",
	"ERROR_IO_FILE_OPEN_CLOSE",
	"ERROR_TREHAD_CREATE",
	"ERROR_THREAD_JOIN",
	"ERROR_META_ANALYSIS",
	"ERROR_PVALUE_FORMAT_ERROR",
	"ERROR_UNDEFINED"
}


#define DONE_NORMAL				(1)
#define DONE_ABNORMAL			(-1)

#ifdef _MSC_VER
	#define	M_PI	(acos(-1.0))	// Accurate PI constant
#endif

// HanEskinPvalueTable Constants
#define TABLE_NROW 			(49)
#define TABLE_MAX_NSTUDY 	(50)
#define TABLE_NCOLUMN 		(331)


class ThreadResult
{
private:
	unsigned int key;
	std::string values_str;

public:
	ThreadResult(unsigned int key_, std::string values_str_)
	{
		this.key = key_;
		this.values_str = values_str_;
	}

	void setKey(unsigned int key_){ this.key = key_};
	void setValues_str(std::string values_str_){ this.values_str = values_str_; }
	unsigned int getKey(void){	return this.key; }
	std::string getValues_str(void){ return this.values_str; }
}


class MetaSnp
{
private:
	// Class-variable part for
	// P-value Table of Han Eskin Statistic.
	// Currently, table dimension is fixed as follows.
	// nStudy (rows) from 2 to 50 (toal 49 rows)
	// statistic thresholds (columns) from 0.0 to 33.0 (total 331 columns)
	static double TABLE_MAX_THRESHOLD;
	static double pvalueTable_[TABLE_NROW][TABLE_NCOLUMN];
	static bool isPvalueTableRead_; // TABLE_NROW x TABLE_NCOLUMN
public:
	static void readPvalueTableFile(std::string pvalueTableFile);

}