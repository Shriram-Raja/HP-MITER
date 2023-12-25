// MPSchedulingPA.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "MPSchedulingFR.h"
#include "RandomSystemsTest.h"
int main(int argc, char * argv[])
{	
	//DebugFR(); return 0;
	MarshallUSweepResult(argc, argv); return 0;
	TestUSweep(argc, argv); return 0;
	//TestOptimizedUnschedCoreComputer(); return 0;		
	//GenerateUSweepTestCases(); return 0;
	//TestUSweep({ 0.5, 4.0 }, 11, "M4N20"); return 0;
	TestMPSchedulingFR();
    return 0;
}

