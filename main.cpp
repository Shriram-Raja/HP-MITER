// MPSchedulingPA.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "MPSchedulingFR.h"
#include "RandomSystemsTest.h"

// Not used in the final experiment
void Example() {
	MPScheduling::RandomSystemParameter param;
	param.cUtilRange = {2.7, 2.7};
	param.cNumProcessorsRange = {4, 4};
	param.cNumTaskRange = { 20, 20 };

	MPScheduling::MultiProcessorSystem cSystem = GenerateRandomSystem(param);

	PriorityAssignment cDkCPA; cDkCPA.GenerateDkCPA(cSystem);
	auto cDkCPAResult = ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(cSystem, cDkCPA);
	cout << "DkC Result: " << cDkCPAResult.first << endl;;

	PriorityAssignment cDAPA(cSystem);
	int iDAResult = DeadlineAnalysis().Audsley(cSystem, cDAPA);
	cout << "DA Result: " << iDAResult << endl;

	MPScheduling::MPSchedulingFR cFR(cSystem);
	cFR.getConfiguration().m_iUseAsyncTest = 0;
	cFR.setCorePerIter(5);
	int iStatus = cFR.FocusRefinement(1, 180);
	auto cStatistic = cFR.GenerateStatisticFile("");
	cStatistic.DisplayStatistic();
}

void Experiment(double util_lb, double util_ub, int num_processors, int num_tasks, int num_point, bool constrained_deadline, const char* folder) {
	// Generate test cases
	auto param = RandomSystemParameter();	
	param.cNumProcessorsRange = {num_processors, num_processors};
	param.cNumTaskRange = {num_tasks, num_tasks};
	param.cPeriodRange = {10, 1000};
	param.iConstrainedDeadline = constrained_deadline;
	const pair<double, double> kUSweepRange = {util_lb, util_ub};	
	
	// Experiment 1
	GenerateUSweepTestCases(param, kUSweepRange, num_point, folder);

	// Experiment 2 for comparison with ML for RT
	// GenerateUSweepTestCases(param, kUSweepRange, num_point, folder, 1);

	TestUSweep(kUSweepRange, num_point, folder);	
	

	// ********************* Aug.5 Add MarshallUSweepResult *********************

	const char* axOutput = "HP_MUTER_8p32_Implicit.txt";
	MarshallUSweepResult(kUSweepRange, num_point, folder, axOutput, 0); 
	
	// ********************* Aug.5 Add MarshallUSweepResult *********************
}

int main(int argc, char * argv[]) {	
	
	{	
	// ********************* Experiment 1 *********************
	// Use MarshallUSweep to generate figure
	// Uncomment the below to enable. n = 4m
	const int kNumPoints = 30;

	// ------------ Implicit Deadline ------------
	// Experiment(0, 4, 4, 16, kNumPoints, false, "HP_MUTER_4p16_Implicit");
	Experiment(0, 8, 8, 32, kNumPoints, false, "HP_MUTER_8p32_Implicit");
	// Experiment(0, 16, 16, 64, kNumPoints, false, "HP_MUTER_16p64_Implicit");
	
	// ------------ Constrained Deadline ------------
	// Experiment(0, 4, 4, 16, kNumPoints, true, "HP_MUTER_4p16_Constrained");
	// Experiment(0, 8, 8, 32, kNumPoints, true, "HP_MUTER_8p32_Constrained");
	// Experiment(0, 16, 16, 64, kNumPoints, true, "HP_MUTER_16p64_Constrained");
	}

	{
	// ********************* Experiment 2 *********************
	// ------------ m = 2, 6 <= n <= 15 ------------
	// Experiment(2, 2, 2, 6, 1, false, "HP_MUTER_Bar_Graph_2p6");
	// Experiment(2, 2, 2, 7, 1, false, "HP_MUTER_Bar_Graph_2p7");
	// Experiment(2, 2, 2, 8, 1, false, "HP_MUTER_Bar_Graph_2p8");
	// Experiment(2, 2, 2, 9, 1, false, "HP_MUTER_Bar_Graph_2p9");
	// Experiment(2, 2, 2, 10, 1, false, "HP_MUTER_Bar_Graph_2p10");
	// Experiment(2, 2, 2, 11, 1, false, "HP_MUTER_Bar_Graph_2p11");
	// Experiment(2, 2, 2, 12, 1, false, "HP_MUTER_Bar_Graph_2p12");
	// Experiment(2, 2, 2, 13, 1, false, "HP_MUTER_Bar_Graph_2p13");
	// Experiment(2, 2, 2, 14, 1, false, "HP_MUTER_Bar_Graph_2p14");
	// Experiment(2, 2, 2, 15, 1, false, "HP_MUTER_Bar_Graph_2p15");

	// ------------ m = 4, 11 <= n <= 20 ------------
	// Experiment(4, 4, 4, 11, 1, false, "HP_MUTER_Bar_Graph_4p11");
	// Experiment(4, 4, 4, 12, 1, false, "HP_MUTER_Bar_Graph_4p12");
	// Experiment(4, 4, 4, 13, 1, false, "HP_MUTER_Bar_Graph_4p13");
	// Experiment(4, 4, 4, 14, 1, false, "HP_MUTER_Bar_Graph_4p14");
	// Experiment(4, 4, 4, 15, 1, false, "HP_MUTER_Bar_Graph_4p15");
	// Experiment(4, 4, 4, 16, 1, false, "HP_MUTER_Bar_Graph_4p16");
	// Experiment(4, 4, 4, 17, 1, false, "HP_MUTER_Bar_Graph_4p17");
	// Experiment(4, 4, 4, 18, 1, false, "HP_MUTER_Bar_Graph_4p18");
	// Experiment(4, 4, 4, 19, 1, false, "HP_MUTER_Bar_Graph_4p19");
	// Experiment(4, 4, 4, 20, 1, false, "HP_MUTER_Bar_Graph_4p20");

	// ------------ m = 6, 16 <= n <= 25 ------------
	// Experiment(6, 6, 6, 16, 1, false, "HP_MUTER_Bar_Graph_6p16");
	// Experiment(6, 6, 6, 17, 1, false, "HP_MUTER_Bar_Graph_6p17");
	// Experiment(6, 6, 6, 18, 1, false, "HP_MUTER_Bar_Graph_6p18");
	// Experiment(6, 6, 6, 19, 1, false, "HP_MUTER_Bar_Graph_6p19");
	// Experiment(6, 6, 6, 20, 1, false, "HP_MUTER_Bar_Graph_6p20");
	// Experiment(6, 6, 6, 21, 1, false, "HP_MUTER_Bar_Graph_6p21");
	// Experiment(6, 6, 6, 22, 1, false, "HP_MUTER_Bar_Graph_6p22");
	// Experiment(6, 6, 6, 23, 1, false, "HP_MUTER_Bar_Graph_6p23");
	// Experiment(6, 6, 6, 24, 1, false, "HP_MUTER_Bar_Graph_6p24");
	// Experiment(6, 6, 6, 25, 1, false, "HP_MUTER_Bar_Graph_6p25");
	}

    return 0;
}

