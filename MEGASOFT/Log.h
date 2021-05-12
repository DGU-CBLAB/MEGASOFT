#pragma once
class Log
{
protected:
	std::string m_argsSummary;
	virtual VOID printLog(DOUBLE time) = 0;
public:
	
};