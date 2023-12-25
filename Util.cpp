#include "stdafx.h"
#include "Util.h"
#include <deque>
#include <fstream>
#include <time.h>
#include "stdafx.h"
#include "Util.h"
#include <assert.h>
#include <pthread.h>
#include <unistd.h>


Timer::Timer()
{
}


Timer::~Timer()
{
}

void Timer::Start()
{
	int iStatus = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &m_tagStartCPUTime);
	assert(iStatus == 0);
}

void Timer::StartThread()
{
	clockid_t cid;
	int iStatus = pthread_getcpuclockid(pthread_self(), &cid);
	assert(iStatus == 0);
	iStatus = clock_gettime(cid, &m_tagStartThreadCPUTime);
        assert(iStatus == 0);
}

void Timer::Stop()
{
	int iStatus = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &m_tagEndCPUTime);
        assert(iStatus == 0);
}

void Timer::StopThread()
{
	clockid_t cid;
        int iStatus = pthread_getcpuclockid(pthread_self(), &cid);
        assert(iStatus == 0);
        iStatus = clock_gettime(cid, &m_tagEndThreadCPUTime);
        assert(iStatus == 0);
}

double Timer::getElapsedTime_s()
{
	double dSec = m_tagEndCPUTime.tv_sec - m_tagStartCPUTime.tv_sec +  (double)(m_tagEndCPUTime.tv_nsec - m_tagStartCPUTime.tv_nsec)/(double)1e9;
	return dSec;
}

double Timer::getElapsedTime_ms()
{
	return 1000 * getElapsedTime_s();
}

double Timer::getElapsedThreadTime_s()
{
        double dSec = m_tagEndThreadCPUTime.tv_sec - m_tagStartThreadCPUTime.tv_sec + (double)(m_tagEndThreadCPUTime.tv_nsec - m_tagStartThreadCPUTime.tv_nsec)/(double)1e9;
        return dSec;
}

double Timer::getElapsedThreadTime_ms()
{
        return getElapsedThreadTime_s() * 1e3;
}

void my_assert(bool bCondition, const char axBuffer[])
{
	if (!bCondition)
	{
		cout << axBuffer << endl;
		while (1);
	}
}

bool IsFileExists(const char axBuffer[])
{
	ifstream ifstreamTest(axBuffer);
	if (ifstreamTest.is_open())
	{
		ifstreamTest.close();
		return true;
	}
	else
	{
		return false;
	}
}

random_device rd;
default_random_engine eng(rd());

double getRandomRealValue(double lb, double ub)
{
	uniform_real_distribution<double> distribution(lb, ub);
	return distribution(eng);
}

int getRandomIntValue(int lb, int ub)
{
	uniform_int_distribution<int> distribution(lb, ub);
	return distribution(eng);
}

bool ProbabilityHit(double dProb)
{
	return getRandomRealValue(0.0, 1.0) <= dProb;
}

double roundReal(double dValue, double iGranularity)
{
	dValue = round(dValue / iGranularity) * iGranularity;
	return dValue;
}

void Makepath(const char axPath[])
{
	char axBuffer[2048] = { 0 };
	sprintf(axBuffer, "mkdir -p %s", axPath);
	system(axBuffer);
}

void Makepath(const string cPath)
{
	Makepath(cPath.data());
}

long long filesize(const char axFileName[]){
	char * axBuffer = new char[1 << 25];
	long long iSize = 0;
	ifstream cFile(axFileName, ios::binary);
	while (cFile.read(axBuffer, 1 << 20) || cFile.gcount() != 0)
		iSize += cFile.gcount();
	cFile.close();
	delete axBuffer;
	return iSize;
};

int ExtractArgument(int argc, char * argv[], unordered_map<string, string> & rmapName2Arg)
{
	char axBuffer[256] = { 0 };
	int iDefaultIndex = 0;
	for (int i = 0; i < argc; i++)
	{
		string stringCommmand = string(argv[i]);
		if (argv[i][0] == '-')
		{
			//name specified
			int iEqualityIndex = stringCommmand.find('=', 1);
			if (iEqualityIndex == -1)
			{
				cout << "I don't recognize this: " << stringCommmand.data() << endl;
				my_assert(0, "");
			}
			string stringArgName = stringCommmand.substr(1, iEqualityIndex - 2 + 1);
			string stringArg = stringCommmand.substr(iEqualityIndex + 1, stringCommmand.size() - iEqualityIndex - 1 + 1);
			rmapName2Arg[stringArgName] = stringArg;
		}
		else
		{
			sprintf(axBuffer, "arg%d", iDefaultIndex);
			iDefaultIndex++;
			string stringArgName(axBuffer);
			rmapName2Arg[stringArgName] = stringCommmand;
		}
	}
	return 0;
}

string itos(long long val) {
	char axBuffer[32] = { 0 };
	sprintf(axBuffer, "lld", val);
	return string(axBuffer);
}

unordered_map<string, string> ExtractArgument(int argc, char * argv[])
{
	unordered_map<string, string> rmapName2Arg;
	char axBuffer[256] = { 0 };
	int iDefaultIndex = 0;
	for (int i = 0; i < argc; i++)
	{
		string stringCommmand = string(argv[i]);
		if (argv[i][0] == '-')
		{
			//name specified
			int iEqualityIndex = stringCommmand.find('=', 1);
			if (iEqualityIndex == -1)
			{
				cout << "I don't recognize this: " << stringCommmand.data() << endl;
				my_assert(0, "");
			}
			string stringArgName = stringCommmand.substr(1, iEqualityIndex - 2 + 1);
			string stringArg = stringCommmand.substr(iEqualityIndex + 1, stringCommmand.size() - iEqualityIndex - 1 + 1);
			rmapName2Arg[stringArgName] = stringArg;
		}
		else
		{
			sprintf(axBuffer, "arg%d", iDefaultIndex);
			iDefaultIndex++;
			string stringArgName(axBuffer);
			rmapName2Arg[stringArgName] = stringCommmand;
		}
	}
	return rmapName2Arg;
}

double getWallTimeSecond()
{
	return (double)clock() / (double)CLOCKS_PER_SEC;
}

string itoa(int iNum) {
	char axBuffer[128];
	sprintf(axBuffer, "%d", iNum);
	return string(axBuffer);
}
string ftoa(double dValue) {
	char axBuffer[128];
	sprintf(axBuffer, "%f", dValue);
	return string(axBuffer);
}

void Sleep(unsigned int millis) {
	usleep(millis * 1000);
}