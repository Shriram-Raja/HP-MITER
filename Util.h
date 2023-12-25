#pragma once
#include <iostream>
using namespace std;
#include <string>
#include <deque>
#include <random>
#include <unordered_map>
#include <time.h>

#define MY_MAX(a, b) (((a) < (b)) ? (b) : (a))
#define MY_MIN(a, b) (((a) > (b)) ? (b) : (a))
#define CLIP(a, b, c) (MY_MIN(MY_MAX(a, b), c))

class StringCompObj
{
public:
	bool operator () (const string & stringLhs, const string & stringRhs) const
	{
		return std::lexicographical_compare(stringLhs.begin(), stringLhs.end(), stringRhs.begin(), stringRhs.end());
	}
};

class DoubleCompObj
{
protected:
	double m_dTolerance;
public:
	DoubleCompObj()
	{
		m_dTolerance = 1e-9;
	}

	DoubleCompObj(double dTol)
	{
		m_dTolerance = dTol < 0 ? -dTol : dTol;
	}

	~DoubleCompObj() {}

	bool operator () (const double & dLHS, const double & dRHS) const
	{
		return isLT(dLHS, dRHS);
	}

	bool isLT(const double & dLHS, const double & dRHS) const
	{
		return (dLHS - dRHS < -m_dTolerance);
	}

	bool isGT(const double & dLHS, const double & dRHS) const
	{
		return (dLHS - dRHS > m_dTolerance);
	}

	bool isGE(const double & dLHS, const double & dRHS) const
	{
		return (dLHS - dRHS >= -m_dTolerance);
	}

	bool isLE(const double & dLHS, const double & dRHS) const
	{
		return (dLHS - dRHS <= m_dTolerance);
	}

	bool isEQ(const double & dLHS, const double & dRHS) const
	{
		return ((dLHS - dRHS) <= m_dTolerance) && ((dLHS - dRHS) >= -m_dTolerance);
	}
};

#define DOUBLE_EQUAL(a,b) (DoubleCompObj().isEQ(a, b))
#define DOUBLE_EQUAL(a,b, tol) (DoubleCompObj(tol).isEQ(a, b))

void my_assert(bool bCondition, const char axBuffer[]);

bool IsFileExists(const char axBuffer[]);

double getRandomRealValue(double lb, double ub);
int getRandomIntValue(int lb, int ub);
bool ProbabilityHit(double dProb);
double roundReal(double dValue, double iGranularity);
void Makepath(const char axPath[]);
void Makepath(const string cPath);

class Timer
{
	struct timespec m_tagStartCPUTime;
	struct timespec m_tagEndCPUTime;
	struct timespec m_tagStartThreadCPUTime;
	struct timespec m_tagEndThreadCPUTime;
public:
	Timer();
	~Timer();
	void Start();
	void Stop();
	void StartThread();
	void StopThread();
	double getElapsedTime_s();
	double getElapsedTime_ms();
	double getElapsedThreadTime_s();
	double getElapsedThreadTime_ms();
	unsigned long long getElapsedTime_100ns();
	unsigned long long getElapsedThreadTime_100ns();
};
long long filesize(const char axFileName[]);
int ExtractArgument(int argc, char * argv[], unordered_map<string, string> & rmapName2Arg);
unordered_map<string, string> ExtractArgument(int argc, char * argv[]);
double getWallTimeSecond();

string itoa(int iNum);
string ftoa(double dValue);

void Sleep(unsigned int millis);
