#include "stdafx.h"
#include "RandomSystemsTest.h"
#include "MPSchedulingFR.h"
#include <fstream>
#include <omp.h>
#include <chrono>
#include <cfloat>

#define CASES 100
#define TIMEOUT 600

using namespace std;
using namespace std::chrono;

namespace MPScheduling {
	RandomSystemParameter::RandomSystemParameter() {
		cUtilRange = { 0.5, 0.5 };
		cNumProcessorsRange = { 4, 4 };
		cNumTaskRange = { 10, 10 };
		cPeriodRange = { 100, 100000 };
		iConstrainedDeadline = 0;
	}

	MultiProcessorSystem GenerateRandomSystem(RandomSystemParameter & rcParam, int MLforRT) {
		double dUtil = getRandomRealValue(rcParam.cUtilRange.first, rcParam.cUtilRange.second);
		int iNumProcessor = getRandomIntValue(rcParam.cNumProcessorsRange.first, rcParam.cNumProcessorsRange.second);
		int iNumTask = getRandomIntValue(rcParam.cNumTaskRange.first, rcParam.cNumTaskRange.second);
		if(MLforRT == 0)
		{
			if(rcParam.iConstrainedDeadline == 0)
				return RandomSystemGenerator().GenerateSystem(iNumTask, iNumProcessor, dUtil, rcParam.cPeriodRange.first, rcParam.cPeriodRange.second);
			else
				return RandomSystemGenerator().GenerateSystemConstrainedDeadline(iNumTask, iNumProcessor, dUtil, rcParam.cPeriodRange.first, rcParam.cPeriodRange.second);
		}
		else
		{
			// Generates random system as described in MLforRT paper
			return RandomSystemGenerator().GenerateSystemMLforRT(iNumTask, iNumProcessor, dUtil, rcParam.cPeriodRange.first, rcParam.cPeriodRange.second);
		}
	}

	void GenerateUSweepTestCasesMT(RandomSystemParameter & rcParam, pair<double, double> cUSweepRange, int iNumPoint, const char axFolder[], int iThreadNum) {		
		vector<thread> vecThreads(iThreadNum);
		double dStep = (cUSweepRange.second - cUSweepRange.first) / (iNumPoint - 1), u = cUSweepRange.first;
		int iCount = 0;
		mutex cULock, cPrintLock;
		for (int i = 0; i < iThreadNum; i++) {
			vecThreads[i] = thread([&]()->void {	
				while (true) {
					cULock.lock();
					double dThisU = u; u += dStep;
					cULock.unlock();
					RandomSystemParameter cCopyParam = rcParam;
					cCopyParam.cUtilRange = { dThisU, dThisU };
					string stringUFolder = string(axFolder) + "/Util" + ftoa(dThisU);
					for (int i = 0; i < 1000; i++) {
						string stringCaseFolder = stringUFolder + "/Case" + itoa(i);
						string stringFileName = stringCaseFolder + "/System.txt";
						if (IsFileExists(stringFileName.data())) continue;
						Makepath(stringCaseFolder.data());
						auto cSystem = GenerateRandomSystem(rcParam);
						cSystem.Write(stringFileName.data());
						cPrintLock.lock();
						printf("\rProgress %d / %d        ", iCount++, iNumPoint * 1000);
						cPrintLock.unlock();
					}
				}				
			});
		}		
		cout << endl;
		for (auto & ele: vecThreads) ele.join();
	}

	void GenerateUSweepTestCases(RandomSystemParameter & rcParam, pair<double, double> cUSweepRange, int iNumPoint, const char axFolder[], int MLforRT = 0) {
		char axBuffer[2048] = { 0 };
		double dStep = 0;
		if(iNumPoint == 1) // Test is performed only on one utilization
		{
			dStep = 1;	
		}
		else
		{
			dStep = (cUSweepRange.second - cUSweepRange.first) / (iNumPoint - 1);
		}
	
		for (double u = cUSweepRange.first; u <= cUSweepRange.second; u += dStep) {
			rcParam.cUtilRange = { u, u };
			string stringUFolder = string(axFolder) + "/Util" + ftoa(u);
			for (int i = 0; i < CASES; i++) {
				printf("\rGenerating u = %.2f, case %d      ", u, i);				
				string stringCaseFolder = stringUFolder + "/Case" + itoa(i);
				string stringFileName = stringCaseFolder + "/System.txt";
				if (IsFileExists(stringFileName.data())) {
					cout << "exist";
					continue;
				}

				cout << "created" << endl;
				Makepath(stringCaseFolder.data());				
				auto cSystem = GenerateRandomSystem(rcParam, MLforRT);
				cSystem.Write(stringFileName.data());
			}
		}
		cout << endl;
	}
	
#ifdef MILP
	void TestUSweep(pair<double, double> cUSweepRange, int iNumPoint, const char axFolder[]) {
		char axBuffer[2048] = { 0 };
		double dStep = 0;
		if(iNumPoint == 1) // Test is performed only on one utilization
		{
			dStep = 1;	
		}
		else
		{
			dStep = (cUSweepRange.second - cUSweepRange.first) / (iNumPoint - 1);
		}
		
		#pragma omp parallel for collapse(2)
		for (int u_count = 0; u_count < iNumPoint; ++u_count) {
		// for (double u = cUSweepRange.first; u <= cUSweepRange.second; u += dStep) {			
			for (int i = 0; i < CASES; i++) {
				double u = cUSweepRange.first + (dStep * (double)(u_count));
				string stringUFolder = string(axFolder) + "/Util" + ftoa(u);
				string stringCaseFolder = stringUFolder + "/Case" + itoa(i);				
				string stringFileName = stringCaseFolder + "/System.txt";
				string stringResultFileName = stringCaseFolder + "/Result.txt";
				if (IsFileExists(stringResultFileName.data())) {
					printf("\rSkip u = %.2f, case %d\n", u, i);
					continue;
				}				
				MultiProcessorSystem cSystem; cSystem.Read(stringFileName.data());

				cout << "----------------------------------" << endl;
				printf("\rTesting u = %.2f, case %d\n", u, i);
				cout << "Case " << i << endl;

				//--------------------------- 1. DkC Heuristic ---------------------------
				//DkC + EPE2
				// Start Timer
				double start0 = getWallTimeSecond(); //omp_get_wtime();
				PriorityAssignment cDkCPA; cDkCPA.GenerateDkCPA(cSystem);
				auto cEPE2PAResult = ResponseTimeAnalysis_EPE().ComputeRTAGSYY2(cSystem, cDkCPA);
				// End Timer
				double end0 = getWallTimeSecond(); //omp_get_wtime();
				double DkC_rt = end0 - start0;
				cout << "DkC + EPE2 Result: " << cEPE2PAResult.first << endl;

				//DkC + GSYY
				auto cDkCPAResult = ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(cSystem, cDkCPA);
				cout << "DkC + GSYY Result: " << cDkCPAResult.first << endl;;	
				
				//DkC + ZLL
				auto cZLLPAResult = ResponseTimeAnalysis_ZLL().ComputeAllRTGivenPA(cSystem, cDkCPA);
				cout << "DkC + ZLL Result: " << cZLLPAResult.first << endl;


				//--------------------------- 2. Deadline Monotonic Priority Ordering ---------------------------
				//DMPO + EPE2
				// Start Timer
				double start1 = getWallTimeSecond(); //omp_get_wtime();
				PriorityAssignment cDMPA; cDMPA.GenerateDMPA(cSystem);
				auto cEPE2DMPAResult = ResponseTimeAnalysis_EPE().ComputeRTAGSYY2(cSystem, cDMPA);
				// End Timer
				double end1 = getWallTimeSecond(); //omp_get_wtime();
				double DMPO_rt = (end1 - start1);
				cout << "DMPO + EPE2 Result: " << cEPE2DMPAResult.first << endl;
				
				//DMPO + GSYY
				auto cDMPAResult = ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(cSystem, cDMPA);
				cout << "DMPO + GSYY Result: " << cDMPAResult.first << endl;
				
				//DMPO + ZLL
				auto cZLLDMPAResult = ResponseTimeAnalysis_ZLL().ComputeAllRTGivenPA(cSystem, cDMPA);
				cout << "ZLL + DMPO Result: " << cZLLDMPAResult.first << endl;


				//--------------------------- 3. Deadline minus Computation time Monotonic Priority Ordering ---------------------------
				//DCMPO + EPE2
				// Start Timer
				double start2 = getWallTimeSecond(); //omp_get_wtime();
				PriorityAssignment cDCMPA; cDCMPA.GenerateDCMPA(cSystem);
				auto cEPE2DCMPAResult = ResponseTimeAnalysis_EPE().ComputeRTAGSYY2(cSystem, cDCMPA);
				// End Timer
				double end2 = getWallTimeSecond(); //omp_get_wtime();
				double DCMPO_rt = (end2 - start2);
				cout << "DCMPO + EPE2 Result: " << cEPE2DCMPAResult.first << endl;
				
				// DCMPO + GSYY
				auto cDCMPAResult = ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(cSystem, cDCMPA);
				cout << "DCMPO + GSYY Result: " << cDCMPAResult.first << endl;
				
				//DCMPO + ZLL
				auto cZLLDCMPAResult = ResponseTimeAnalysis_ZLL().ComputeAllRTGivenPA(cSystem, cDCMPA);
				cout << "ZLL + DCMPO Result: " << cZLLDCMPAResult.first << endl;
				
				//--------------------------- 4. Audsley's Optimal Priority Assignment ---------------------------
				// Start Timer
				double start3 = getWallTimeSecond(); //omp_get_wtime();
				PriorityAssignment cDAPA(cSystem);
				int iDAResult = DeadlineAnalysis().Audsley(cSystem, cDAPA, DA);
				// End Timer
				double end3 = getWallTimeSecond(); //omp_get_wtime();
				double DA_rt = (end3 - start3);
				cout << "DA Result: " << iDAResult << endl;

				//--------------------------- 5. MUTER ---------------------------
				MPSchedulingFR cFR(cSystem);
				cFR.setCorePerIter(5);
				// Start Timer
				double MUTER_start = getWallTimeSecond(); //omp_get_wtime();
				int iStatus = cFR.FocusRefinement(1, TIMEOUT);
				// End Timer
				double MUTER_time = getWallTimeSecond() - MUTER_start;
				double MUTER_rt = MUTER_time;
				auto cStatistic = cFR.GenerateStatisticFile("");
				cout << "FR Result: " << iStatus << endl;
				double MUTER_timeout = (iStatus == -1) ? MUTER_time : 0;
				double MUTER_unsched = (iStatus == 0) ? MUTER_time : 0;
				double MUTER_sched = (iStatus == 1) ? MUTER_time : 0;
				assert((iDAResult != 1) || (iStatus == 1));
				cout << "MUTER: " << MUTER_time 
					 << "; MUTER_timeout: " << MUTER_timeout
					 << "; MUTER_unsched: " << MUTER_unsched
					 << "; MUTER_sched: " << MUTER_sched 
				     << "; MUTER Stat file: " << cStatistic.getItem("Total Time (ms)").getValue() << endl;

				//--------------------------- 6. Hybrid Priority Assignment Algorithm with MUTER ---------------------------
				PriorityAssignment cHP_MUTER;
				// Start Timer
				double start5 = getWallTimeSecond(); //omp_get_wtime();
				auto HP_MUTER_Result = cHP_MUTER.GenerateHP_MUTER(cSystem, 0.40, false, true);
				// End Timer
				double end5 = getWallTimeSecond(); //omp_get_wtime();
				double HP_MUTER_rt = (end5 - start5);
				cout << "New PA Result: " << HP_MUTER_Result << endl;
				double HP_MUTER_unsched = (HP_MUTER_Result == 0) ? HP_MUTER_rt : 0;
				double HP_MUTER_sched = (HP_MUTER_Result == 1) ? HP_MUTER_rt : 0;
				double HP_MUTER_EarlyExit = (HP_MUTER_Result == 2) ? HP_MUTER_rt : 0;
				cout << "; HP_MUTER: " << HP_MUTER_rt 
					 << "; HP_MUTER_unsched: " << HP_MUTER_unsched
					 << "; HP_MUTER_sched: " << HP_MUTER_sched 
					 << "; HP_MUTER_EarlyExit: " << HP_MUTER_EarlyExit << endl;

				//--------------------------- 7. DA-LC + Audsley's Optimal Priority Assignment ---------------------------
				// Start Timer
				double start6 = getWallTimeSecond(); //omp_get_wtime();
				PriorityAssignment cDALCPA(cSystem);
				int iDA_LCResult = DeadlineAnalysis().Audsley(cSystem, cDALCPA, DA_LC);
				// End Timer
				double end6 = getWallTimeSecond(); //omp_get_wtime();
				double DA_LC_rt = (end6 - start6);
				cout << "DA-LC Result: " << iDA_LCResult << endl;

				//--------------------------- 8. C-RTA + Audsley's Optimal Priority Assignment ---------------------------
				PriorityAssignment cCRTAPA(cSystem);
				int iC_RTAResult = DeadlineAnalysis().Audsley(cSystem, cCRTAPA, C_RTA);
				cout << "C-RTA Result: " << iC_RTAResult << endl;
				
				//--------------------------- 9., 10. HP-MUTER modifications ---------------------------
				PriorityAssignment cHP_MUTER_late; // Step 23, i.e., Lateness removed 
				// Start Timer
				double start7 = getWallTimeSecond(); //omp_get_wtime();
				auto HP_MUTER_Result_late = cHP_MUTER_late.GenerateHP_MUTER(cSystem, 0.40, false, false);
				// End Timer
				double end7 = getWallTimeSecond(); //omp_get_wtime();
				double HP_MUTER_late_rt = (end7 - start7);
				
				PriorityAssignment cHP_MUTER_DkC; // use DkC if MUTER returns infeasible or times out
				// Start Timer
				double start8 = getWallTimeSecond(); //omp_get_wtime();
				auto HP_MUTER_Result_DkC = cHP_MUTER_DkC.GenerateHP_MUTER(cSystem, 0.40, true, true);
				// End Timer
				double end8 = getWallTimeSecond(); //omp_get_wtime();
				double HP_MUTER_DkC_rt = (end8 - start8);
				
				//--------------------------- Write results into the file ---------------------------
				ofstream cOutputFile(stringResultFileName.data());
				
				// 1. DkCPA Output - 0, 1, 2
				cOutputFile << cDkCPAResult.first << endl;
				cOutputFile << cZLLPAResult.first << endl;
				cOutputFile << cEPE2PAResult.first << endl;
				
				// 2. DMPA Output - 3, 4, 5
				cOutputFile << cDMPAResult.first << endl;
				cOutputFile << cZLLDMPAResult.first << endl;
				cOutputFile << cEPE2DMPAResult.first << endl;

				// 3. DCMPA Output - 6, 7, 8
				cOutputFile << cDCMPAResult.first << endl;
				cOutputFile << cZLLDCMPAResult.first  << endl;
				cOutputFile << cEPE2DCMPAResult.first << endl;
				
				// 4. DA - 9
				cOutputFile << iDAResult << endl;

				// 5. MUTER - 10
				cOutputFile << iStatus << endl;

				// 6. HP-MUTER - 11
				cOutputFile << HP_MUTER_Result << endl;

				// MUTER time from Statistic file - 12
				cOutputFile << cStatistic.getItem("Total Time (ms)").getValue() << endl;

				// MUTER time from here - 13, 14, 15, 16
				cOutputFile << MUTER_time << endl;
				cOutputFile << MUTER_timeout << endl;
				cOutputFile << MUTER_unsched << endl;
				cOutputFile << MUTER_sched << endl;

				// HP-MUTER time from here - 17, 18, 19, 20
				cOutputFile << HP_MUTER_rt << endl;
				cOutputFile << HP_MUTER_unsched << endl;
				cOutputFile << HP_MUTER_sched << endl;
				cOutputFile << HP_MUTER_EarlyExit << endl;

				// MLforRT type info - 21, 22
				bool schedAtLeastOne = (cDkCPAResult.first	|| cZLLPAResult.first		|| cEPE2PAResult.first 
									 || cDMPAResult.first	|| cZLLDMPAResult.first		|| cEPE2DMPAResult.first
									 || cDCMPAResult.first	|| cZLLDCMPAResult.first	|| cEPE2DCMPAResult.first
									 || iDAResult 			|| iStatus 					|| HP_MUTER_Result);
				bool schedOnlyHPMUTER = HP_MUTER_Result && !(cDkCPAResult.first	|| cZLLPAResult.first		|| cEPE2PAResult.first 
									 					  || cDMPAResult.first	|| cZLLDMPAResult.first		|| cEPE2DMPAResult.first
									 					  || cDCMPAResult.first	|| cZLLDCMPAResult.first	|| cEPE2DCMPAResult.first
									 					  || iDAResult 			|| iStatus);

				cOutputFile << schedAtLeastOne << endl;
				cOutputFile << schedOnlyHPMUTER << endl;

				// C-RTA + OPA - 23
				cOutputFile << iC_RTAResult << endl;				

				// DA-LC + OPA - 24
				cOutputFile << iDA_LCResult << endl;

				// ZLL Compiled - 25
				auto ZLL_Heu_Compiled = cZLLPAResult.first || cZLLDMPAResult.first || cZLLDCMPAResult.first;
				cOutputFile << ZLL_Heu_Compiled << endl;

				// RTA-LC / GSYY Compiled - 26
				auto RTA_LC_Heu_Compiled = cDkCPAResult.first || cDMPAResult.first || cDCMPAResult.first;
				cOutputFile << RTA_LC_Heu_Compiled << endl;
				
				// EPE Compiled - 27
				auto EPE_Heu_Compiled = cEPE2PAResult.first || cEPE2DMPAResult.first || cEPE2DCMPAResult.first;
				cOutputFile << EPE_Heu_Compiled << endl;

				float task_util_ = 0;
				for(int util_iter = 0; util_iter < cSystem.getNumTasks(); ++util_iter)
				{
					Task T_ui = cSystem.getTask(util_iter);
					task_util_ += ((float)T_ui.C()) / ((float)T_ui.T());
				}

				// obtain non-zero integer, f such that (f - 1)m < util of taskset <= fm
				int factor_util = 1;
				for(; factor_util <= 10; ++factor_util)
				{
					if(task_util_ <= factor_util * 0.1 * cSystem.getNumProcessors())
					{
						break;
					}
				}
				// Utilization Factor - 28
				cOutputFile << factor_util << endl;

				// HP-MUTER mods - 29, 30
				cOutputFile << HP_MUTER_Result_late << endl;
				cOutputFile << HP_MUTER_Result_DkC << endl;
				
				// Times DkC, DMPO, DCMPO, OPA+DA, MUTER, OPA+DA-LC, HP-MUTER mods - 31, 32, 33, 34, 35, 36, 37, 38, 39
				cOutputFile << DkC_rt << endl;
				cOutputFile << DMPO_rt << endl;
				cOutputFile << DCMPO_rt << endl;
				cOutputFile << DA_rt << endl;
				cOutputFile << MUTER_rt << endl;
				cOutputFile << HP_MUTER_rt << endl;
				cOutputFile << DA_LC_rt << endl;
				cOutputFile << HP_MUTER_late_rt << endl;
				cOutputFile << HP_MUTER_DkC_rt << endl;
				
				cOutputFile.close();
			}
		}
	}

	void testSingle(int argc, char * argv[]) {
		MultiProcessorSystem cSystem; cSystem.Read(argv[1]);

		PriorityAssignment cDkCPA; cDkCPA.GenerateDkCPA(cSystem);
		auto cDkCPAResult = ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(cSystem, cDkCPA);
		cout << "DkC Result: " << cDkCPAResult.first << endl;;

		PriorityAssignment cDAPA(cSystem);
		int iDAResult = DeadlineAnalysis().Audsley(cSystem, cDAPA, DA);
		cout << "DA Result: " << iDAResult << endl;

		MPSchedulingFR cFR(cSystem);
		cFR.getConfiguration().m_iUseAsyncTest = 0;
		cFR.setCorePerIter(5);
		int iStatus = cFR.FocusRefinement(0, 180);
		auto cStatistic = cFR.GenerateStatisticFile("");

		system("pause");
	}

	void TestUSweepMT(pair<double, double> cUSweepRange, int iNumPoint, const char axFolder[], int iThreadNum, int iFRNoAsync) {
		double dStep = (cUSweepRange.second - cUSweepRange.first) / (iNumPoint - 1);
		vector<ofstream> vecThreadLogFile(iThreadNum);
		for(int i = 0; i < iThreadNum; i++) vecThreadLogFile[i].open(string(string("Thread") + itoa(i) + "Log.txt").data());
		for (double u = cUSweepRange.first; u <= cUSweepRange.second; u += dStep) {
			cout << "MT testing u = " << u << endl;
			string stringUFolder = string(axFolder) + "/Util" + ftoa(u);
			vector<thread> vecThreads(iThreadNum);			
			volatile int iIndex = 0;
			mutex cLock;
			for (int iThread = 0; iThread < vecThreads.size(); iThread++) {							
				vecThreads[iThread] = thread([&](int iThreadIndex) {
					while (true) {
						cLock.lock();
						int i = iIndex++;
						cLock.unlock();
						if (i >= 1000) break;
						string stringCaseFolder = stringUFolder + "/Case" + itoa(i);
						string stringFileName = stringCaseFolder + "/System.txt";
						string stringResultFileName = iFRNoAsync == 0 ? stringCaseFolder + "/Result.txt" : stringCaseFolder + "/FRNoAsyncResult.txt";
						if (IsFileExists(stringResultFileName.data())) {
							vecThreadLogFile[iThreadIndex] << "Skip u = " << u << ", case " << i << endl;
							continue;
						}
						MultiProcessorSystem cSystem; cSystem.Read(stringFileName.data());

						vecThreadLogFile[iThreadIndex] << "----------------------------------" << endl;
						vecThreadLogFile[iThreadIndex] << "Testing u = " << u << ", case " << i << endl;
						//DkC + GSYY
						PriorityAssignment cDkCPA; cDkCPA.GenerateDkCPA(cSystem);
						auto cDkCPAResult = ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(cSystem, cDkCPA);
						vecThreadLogFile[iThreadIndex] << "DkC Result: " << cDkCPAResult.first << endl;;
						//DA
						PriorityAssignment cDAPA(cSystem);
						int iDAResult = DeadlineAnalysis().Audsley(cSystem, cDAPA, DA);
						vecThreadLogFile[iThreadIndex] << "DA Result: " << iDAResult << endl;
						//FR
						MPSchedulingFR cFR(cSystem);
						cFR.getConfiguration().m_iUseAsyncTest = iFRNoAsync == 0;
						cFR.setCorePerIter(5);
						int iStatus = cFR.FocusRefinement(0, 180);
						auto cStatistic = cFR.GenerateStatisticFile("");
						vecThreadLogFile[iThreadIndex] << "FR Result: " << iStatus << endl;
						if (iDAResult == 1 && iStatus == 0) cout << "FR not optimal? " << endl;
						assert(((iDAResult == 1 && iStatus == 0)) == false);

						ofstream cOutputFile(stringResultFileName.data());
						cOutputFile << cDkCPAResult.first << endl;
						cOutputFile << iDAResult << endl;
						cOutputFile << iStatus << endl;
						cOutputFile << cStatistic.getItem("Total Time (ms)").getValue() << endl;
						cOutputFile.close();
					}
									
				}, iThread);
			}	
			Sleep(500);
			for (auto & ele : vecThreads) ele.join();
		}
		for (auto & ele : vecThreadLogFile) ele.close();
	}	

	void TestUSweep(int argc, char * argv[]) {		
		try {
			auto cArgs = ExtractArgument(argc, argv);
			int iNumberofTasks = atoi(cArgs.at("ntasks").data());
			cout << "ntasks = " << iNumberofTasks << endl;
			int iNumberofProcessors = atoi(cArgs.at("ncpus").data());
			cout << "ncpus = " << iNumberofProcessors << endl;
			double dULB = atof(cArgs.at("dLB").data());
			double dUUB = atof(cArgs.at("dUB").data());
			printf("U: %.2f ~ %.2f\n", dULB, dUUB);
			int iNPoints = atof(cArgs.at("n").data());
			cout << "n = " << iNPoints << endl;
			int iConstrainedDeadline = atoi(cArgs.at("ConstrainedDeadline").data());
			cout << "ConstrainedDeadline = " << iConstrainedDeadline << endl;
			int iNumThreads = atoi(cArgs.at("#Threads").data());https://client.teamviewer.com/taf/creatives/redesign-template2016/img/tv2019-q4-nov-ppeula-en-de-core-0_tv2019-q4-nov-ppeula-global/tv2019-q4-nov-ppeula-en-de-core-0_tv2019-q4-nov-ppeula-global-en.jpg
			cout << "#Threads = " << iNumThreads << endl;
			const char * axFolder = cArgs.at("Folder").data();		
			cout << "Folder = " << axFolder << endl;
			int iFRNoAsynTest = cArgs.count("FRNoAsync") ? atoi(cArgs.at("FRNoAsync").data()) : 0;
			RandomSystemParameter cParam;
			cParam.cNumTaskRange = { iNumberofTasks, iNumberofTasks };
			cParam.cNumProcessorsRange = { iNumberofProcessors, iNumberofProcessors };
			cParam.cPeriodRange = { 100, 100000 };
			cParam.iConstrainedDeadline = iConstrainedDeadline;
			GenerateUSweepTestCases(cParam, { dULB, dUUB }, iNPoints, axFolder);
			TestUSweepMT({ dULB, dUUB }, iNPoints, axFolder, iNumThreads, iFRNoAsynTest);
		}
		catch (std::out_of_range & err) {
			cout << err.what() << endl;
		}		
	}

	
#endif

// Reads all the Result.txt files and compiles the result into a single txt file
	void MarshallUSweepResult(pair<double, double> cUSweepRange, int iNumPoint, const char axFolder[], const char axOutputFileName[], int iFRNoAsync) {
		char axBuffer[2048] = { 0 };
		double dStep = (cUSweepRange.second - cUSweepRange.first) / (iNumPoint - 1);

		ofstream cOutputFile(axOutputFileName);
		assert(cOutputFile.is_open());
		
		for (double u = cUSweepRange.first; u <= cUSweepRange.second; u += dStep) {
			string stringUFolder = string(axFolder) + "/Util" + ftoa(u);
			vector<double> vecResultSum(40,0);
			
			int MUTER_timeout = 0;
			int MUTER_unsched = 0;
			int MUTER_sched = 0;

			int HP_MUTER_unsched = 0;		
			int HP_MUTER_sched = 0;
			int HP_MUTER_EarlyExit = 0;

			int one = 0,	two = 0,
				three = 0,	four = 0, 
				five = 0,	six = 0, 
				seven = 0,	eight = 0, 
				nine = 0, 	ten = 0,
				eleven = 0;

			// max, min, and avg for all priority assignment algorithms
			// 0 DkC, 1 DMPO, 2 DCMPO, 3 DA, 4 MUTER, 5 HP-MUTER, 6 DA-LC, 7 HP-MUTER late, 8 HP-MUTER DkC
			vector<double> vecMax(9, DBL_MIN);
			vector<double> vecMin(9, DBL_MAX);

			bool bValid = true;
			for (int i = 0; i < CASES; i++) {
				string stringCaseFolder = stringUFolder + "/Case" + itoa(i);				
				string stringFileName = stringCaseFolder + "/System.txt";
				MultiProcessorSystem cSystem; cSystem.Read(stringFileName.data());
				//string stringResultFileName = stringCaseFolder + "/Result.txt";
				string stringResultFileName = iFRNoAsync == 0 ? stringCaseFolder + "/Result.txt" : stringCaseFolder + "/FRNoAsyncResult.txt";
				ifstream cInputFile(stringResultFileName.data());
				cout << stringResultFileName.data() << " " << cInputFile.is_open() << endl;;
				bValid = cInputFile.is_open();				
				if (!bValid) break;	// if there isn't enough data break the loop and end 

				vector<double> vecResult(40,0);
				for (int j = 0; j < 40; j++) cInputFile >> vecResult[j];
				
				// Number of schedulable task sets in each PA + SA pair
				for (int j = 0; j < 12; j++) vecResultSum[j] += vecResult[j] > 0;

				// MUTER
				if(vecResult[10] == -1)
				{
					if(vecResult[14] == 0)
					{
						cout << "MUTER_timeout = 0 but result is -1\n";
					}
					MUTER_timeout++;
				}
				if (vecResult[10] == 0)
				{
					if(vecResult[15] == 0)
					{
						cout << "MUTER_unsched = 0 but result is 0\n";
					}
					MUTER_unsched++;
				}
				if (vecResult[10] == 1)
				{
					if(vecResult[16] == 0)
					{
						cout << "MUTER_sched = 0 but result is 1\n";
					}
					MUTER_sched++;
				}

				// HPA with MUTER
				if(vecResult[11] == 0)
				{
					if(vecResult[18] == 0)
					{
						cout << "HP_MUTER_unsched = 0 but result is 0\n";
					}
					HP_MUTER_unsched++;
				}
				if (vecResult[11] == 1)
				{
					if(vecResult[19] == 0)
					{
						cout << "HP_MUTER_sched = 0 but result is 1\n";
					}
					HP_MUTER_sched++;
				}
				if (vecResult[11] == 2)
				{
					if(vecResult[20] == 0)
					{
						cout << "HP_MUTER_EarlyExit = 0 but result is 2\n";
					}
					HP_MUTER_EarlyExit++;
				}
				
				
				// Sum of total time (ms) - from MUTER Statistic file
				vecResultSum[12] += vecResult[12];

				// Sum of total MUTER time from here
				vecResultSum[13] += vecResult[13];
				vecResultSum[14] += vecResult[14];
				vecResultSum[15] += vecResult[15];
				vecResultSum[16] += vecResult[16];

				// Sum of total HP_MUTER time from here
				vecResultSum[17] += vecResult[17];
				vecResultSum[18] += vecResult[18];
				vecResultSum[19] += vecResult[19];
				vecResultSum[20] += vecResult[20];

				for (int z = 21; z < 28; z++) vecResultSum[z] += vecResult[z] > 0;

				// C-RTA Check = !C-RTA && (MUTER or GSYY+Heu Compiled)
				if(!vecResult[23] && (vecResult[10] || vecResult[26]))
				{
					cout << "C-RTA does not dominate RTA-LC" << endl;
				}

				// DA-LC Check = !DA-LC && DA
				if(!vecResult[24] && vecResult[9])
				{
					cout << "DA-LC does not dominate DA" << endl;
				}

				switch ((int)(vecResult[28]))
				{
					case 1:
						++one;
						break;

					case 2:
						++two;
						break;

					case 3:
						++three;
						break;
					
					case 4:
						++four;
						break;

					case 5:
						++five;
						break;

					case 6:
						++six;
						break;
					
					case 7:
						++seven;
						break;

					case 8:
						++eight;
						break;

					case 9:
						++nine;
						break;

					case 10:
						++ten;
						break;
				
					default:
						++eleven;
						cout << "ELEVEN!! \n CASE : " << i << endl;
						break;
				}

				// number of sched HP-MUTER mods late, DkC
				for (int x = 29; x < 31; x++) vecResultSum[x] += vecResult[x] > 0;

				// sum for all priority assignment algorithms
				// 31 DkC, 32 DMPO, 33 DCMPO, 34 DA, 35 MUTER, 36 HP-MUTER, 37 DA-LC, 38 HP-MUTER late, 39 HP-MUTER DkC
				for (int x = 31; x < 40; x++) vecResultSum[x] += vecResult[x];

				cout << "Case " << i << " ";
				if(vecResult[11] > 0)
				{
					cout << "HP_MUTER = " << vecResult[11] << "\n";
				}
				if(vecResult[29] > 0)
				{
					cout << "HP_MUTER_late = " << vecResult[29] << "\n";
				}
				if(vecResult[30] > 0)
				{
					cout << "HP_MUTER_DkC = " << vecResult[30] << "\n";
				}

				// max and min for all priority assignment algorithms
				// 0 DkC, 1 DMPO, 2 DCMPO, 3 DA, 4 MUTER, 5 HP-MUTER, 6 DA-LC, 7 HP-MUTER late, 8 HP-MUTER DkC
				for(int y = 0; y < 9; ++y)
				{
					vecMax[y] = max(vecMax[y], vecResult[y + 31]);
					vecMin[y] = min(vecMin[y], vecResult[y + 31]);
				}

				cInputFile.close();
			}
			if (bValid) {
				cOutputFile << u << " ";

				for (int j = 0; j < 21; j++) cOutputFile << vecResultSum[j] << " ";

				cOutputFile << MUTER_timeout << " ";
				cOutputFile << MUTER_unsched << " ";
				cOutputFile << MUTER_sched << " ";
				cOutputFile << HP_MUTER_unsched << " ";
				cOutputFile << HP_MUTER_sched << " ";
				cOutputFile << HP_MUTER_EarlyExit << " ";
				
				for (int j = 21; j < 28; j++) cOutputFile << vecResultSum[j] << " ";
				
				// cOutputFile << one << " ";
				// cOutputFile << two << " ";
				// cOutputFile << three << " ";
				// cOutputFile << four << " ";
				// cOutputFile << five << " ";
				// cOutputFile << six << " ";
				// cOutputFile << seven << " ";
				// cOutputFile << eight << " ";
				// cOutputFile << nine << " ";
				// cOutputFile << ten << " ";
				// cOutputFile << eleven << " ";

				for (int j = 29; j < 40; j++) cOutputFile << vecResultSum[j] << " ";

				for(int y = 0; y < 9; ++y)
				{
					cOutputFile << vecMax[y] << " ";
					cOutputFile << vecMin[y] << " ";
				}

				cOutputFile << endl;
			}			
		}
		cOutputFile.close();
	}
	
	void MarshallUSweepResult(int argc, char * argv[]) {
		auto cArgs = ExtractArgument(argc, argv);		
		double dULB = atof(cArgs.at("dLB").data());
		double dUUB = atof(cArgs.at("dUB").data());
		int iNPoints = atof(cArgs.at("n").data());
		const char * axFolder = cArgs.at("Folder").data();
		const char * axOutputFile = cArgs.at("Output").data();
		printf("{%.2f, %.2f}, %d, %s, %s\n", dULB, dUUB, iNPoints, axFolder, axOutputFile);
		int iFRNoAsynTest = cArgs.count("FRNoAsync") ? atoi(cArgs.at("FRNoAsync").data()) : 0;
		MarshallUSweepResult({ dULB, dUUB }, iNPoints, axFolder, axOutputFile, iFRNoAsynTest);
	}	
}
