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

#define ERROR_IO				(100)
#define ERROR_IO_FILE_CLOSE		(101)
#define ERROR_TREHAD_CREATE 	(200)
#define ERROR_THREAD_JOIN		(201)
#define ERROR_META_ANALYSIS		(300)


#define DONE_NORMAL				(1)
#define DONE_ABNORMAL			(-1)

#ifdef _MSC_VER
	#define	M_PI	(acos(-1.0))	// Accurate PI constant
#endif

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
