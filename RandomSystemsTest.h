#pragma once
#include "MPSchedulingFR.h"

namespace MPScheduling {
	struct RandomSystemParameter {
		pair<double, double> cUtilRange;
		pair<int, int> cNumProcessorsRange;
		pair<int, int> cNumTaskRange;		
		pair<ValueType, ValueType> cPeriodRange;
		int iConstrainedDeadline;
		RandomSystemParameter();
	};

	MultiProcessorSystem GenerateRandomSystem(RandomSystemParameter & rcParam, int MLforRT = 0);
	void GenerateUSweepTestCases(RandomSystemParameter & rcParam, pair<double, double> cUSweepRange, int iNumPoint, const char axFolder[], int MLforRT = 0);
	void GenerateUSweepTestCasesMT(RandomSystemParameter & rcParam, pair<double, double> cUSweepRange, int iNumPoint, const char axFolder[], int iThreadNum);
	void GenerateUSweepTestCases();
	void TestUSweep(pair<double, double> cUSweepRange, int iNumPoint, const char axFolder[]);
	void TestUSweepMT(pair<double, double> cUSweepRange, int iNumPoint, const char axFolder[], int iThreadNum, int iFRNoAsync);

	// add the following declaration
	void MarshallUSweepResult(pair<double, double> cUSweepRange, int iNumPoint, const char axFolder[], const char axOutputFileName[], int iFRNoAsync);
	//
	void MarshallUSweepResult(int argc, char * argv[]);
	void TestUSweep(int argc, char * argv[]);		
	void testSingle(int argc, char * argv[]);
}
