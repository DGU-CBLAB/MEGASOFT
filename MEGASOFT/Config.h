#pragma once

#ifndef _CONFIG
#define _CONFIG

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include<iostream>
#include<string>

#define VOID	void
#define INT32	int
#define UINT32	unsigned int
#define DOUBLE	double
#define LONG	long
#define BOOL	bool

#define TRUE	true
#define FALSE	false

#define FILE_READ	("r")
#define FILE_WRITE	("w")

#define DELIM_TAB	"\t"
#define DELIM_SPACE	" "

#define CHAR_SHARP '#'

#define SECONDS_IN_MINUTE	(60.0)

#define ERROR_IO				100
#define ERROR_IO_FILE_CLOSE		101
#define ERROR_THREAD_CREATE		200
#define ERROR_THREAD_JOIN		201
#define ERROR_META_ANALYSIS		300
#define ERROR_FILE_NOT_FOUND	404
#define ERROR_FILE_READ_FAIL	405

#define DONE_NORMAL			1
#define DONE_ABNORMAL		-1
#define ERR_THREAD_CREATE	100
#define ERR_THREAD_JOIN		101

/// DEFAULT OPTION VALUES
#define DEFAULT_OUTPUT_FILE					("out")
#define DEFAULT_PVALUE_TABLE_FILE			("HanEskinPvalueTable.txt")
#define DEFAULT_LOG_FILE					(".\\log")
#define DEFAULT_INPUT_LAMBDA_MEAN_EFFECT	(1.0)
#define DEFAULT_INPUT_LAMBDA_HETEROGENEITY	(1.0)
#define DEFAULT_WILL_COMPUTE_MVALUE			(FALSE)
#define DEFAULT_PRIOR_SIGMA					(0.2)
#define DEFAULT_PRIOR_ALPHA					(1.0)
#define DEFAULT_PRIOR_BETA					(1.0)
#define DEFAULT_MVALUE_PVALUE_THRESHOLD		(1E-7)
#define DEFAULT_MVALUE_METHOD				("exact")
#define DEFAULT_MCMC_SAMPLES				(10000)
#define DEFAULT_MCMC_BURNIN					(1000)
#define DEFAULT_MCMC_PROB_RANDOM			(0.01)
#define DEFAULT_MCMC_MAX_NUM_FLIP			(0.1)
#define DEFAULT_WILL_COMPUTE_BINARY_EFFECTS (FALSE)
#define DEFAULT_BINARY_EFFECTS_SAMPLES		(1000)
#define DEFAULT_BINARY_EFFECTS_LARGE_SAMPLE (100000)
#define DEFAULT_BINARY_EFFECTS_PVALUE_THRESHOLD	(1E-4)
#define DEFAULT_SEED						(0)
#define DEFAULT_VERBOSE						(FALSE)
#define DEFAULT_CPU_NUM_THREAD				(1)
#define DEFAULT_GPU							(FALSE)

// HanEskin Table
#define HANESKIN_TABLE_ROW	(49)
#define HANESKIN_TABLE_COL	(331)

#endif

