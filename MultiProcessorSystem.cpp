#include "stdafx.h"
#include "MultiProcessorSystem.h"
#include <assert.h>
#include <random>
#include "Util.h"
#include <fstream>
#include <algorithm>
#include "MPSchedulingFR.h"
#include <limits>

#define DEBUG 0
#define DEBUG_FINE 0
#define DEBUG_PROP 0

#define PROP_CHANGE 0

#define TIMEOUT 600

namespace MPScheduling {
	Task::Task():m_dWCET(0), m_dPeriod(0), m_dDeadline(0) {}

	Task::Task(ValueType dWCET, ValueType dPeriod, ValueType dDeadline)
	:m_dWCET(dWCET), m_dPeriod(dPeriod), m_dDeadline(dDeadline){

	}

	ostream & operator << (ostream & os, const Task & rcTask) {
		os << rcTask.C() << ' ' << rcTask.D() << ' ' << rcTask.T();
		return os;
	}

	istream & operator >> (istream & is, Task & rcTask) {
		ValueType T, D, C;
		is >> C >> D >> T;
		rcTask = Task(C, T, D);
		return is;
	}

	MultiProcessorSystem::MultiProcessorSystem(int iNumProcessors): m_iNumProcessors(iNumProcessors){
	}

	MultiProcessorSystem::~MultiProcessorSystem(){
	}

	void MultiProcessorSystem::Write(ofstream & cOutputFile) {
		cOutputFile << m_iNumProcessors << ' ' << m_vectorTasks.size() << endl;		
		for (auto & ele : m_vectorTasks) cOutputFile << ele << endl;
	}

	void MultiProcessorSystem::Read(ifstream & cInputFile) {
		int iNumTasks = 0;
		cInputFile >> m_iNumProcessors >> iNumTasks;
		m_vectorTasks.reserve(iNumTasks);
		for (int i = 0; i < iNumTasks; i++) {
			m_vectorTasks.push_back(Task());
			cInputFile >> m_vectorTasks.back();
		}
	}

	void MultiProcessorSystem::Read(const char axFileName[]) {
		ifstream cInputFile(axFileName);
		assert(cInputFile.is_open());
		Read(cInputFile);
		cInputFile.close();
	}

	void MultiProcessorSystem::Write(const char axFileName[]) {
		ofstream cOutputFile(axFileName);
		assert(cOutputFile.is_open());
		Write(cOutputFile);
		cOutputFile.close();
	}

	double MultiProcessorSystem::TotalUtil() const {
		double ret = 0;
		for (auto & ele : m_vectorTasks) ret += ele.U();
		return ret;
	}

	PriorityAssignment::PriorityAssignment(MultiProcessorSystem & rcSystem)
	: m_vectorTask2Priority(rcSystem.getNumTasks(), -1), m_vectorPriority2Task(rcSystem.getNumTasks(), -1){
		
	}

	int PriorityAssignment::getPriority(int iTaskIndex) const{
		assert(iTaskIndex < m_vectorTask2Priority.size());
		return m_vectorTask2Priority[iTaskIndex];
	}

	int PriorityAssignment::getTask(int iPriority) const{
		assert(iPriority < m_vectorPriority2Task.size());
		return m_vectorPriority2Task[iPriority];
	}

	int PriorityAssignment::getSize() const {
		int iCount = 0;
		for (auto ele : m_vectorPriority2Task) iCount += ele != -1;
		return iCount;
	}

	void PriorityAssignment::setPriority(int iTaskIndex, int iPriority) {
		assert(getTask(iPriority) == -1);
		assert(getPriority(iTaskIndex) == -1);
		m_vectorPriority2Task[iPriority] = iTaskIndex;
		m_vectorTask2Priority[iTaskIndex] = iPriority;
	}

	void PriorityAssignment::unset(int iTaskIndex) {
		if (iTaskIndex == -1 || getPriority(iTaskIndex) == -1) return;
		m_vectorPriority2Task[getPriority(iTaskIndex)] = -1;
		m_vectorTask2Priority[iTaskIndex] = -1;
	}

	void PriorityAssignment::Write(const char axFileName[]) {
		ofstream cOutputFile(axFileName);
		assert(cOutputFile.is_open());
		cOutputFile << m_vectorTask2Priority.size() << endl;
		for (auto ele : m_vectorTask2Priority) cOutputFile << ele << endl;
		cOutputFile.close();
	}

	void PriorityAssignment::Read(const char axFileName[]) {
		ifstream cInputFile(axFileName);
		assert(cInputFile.is_open());
		int iNumTasks = 0;
		cInputFile >> iNumTasks;
		m_vectorPriority2Task = vector<int>(iNumTasks, -1);
		m_vectorTask2Priority = m_vectorPriority2Task;
		for (int i = 0; i < iNumTasks; i++) {
			cInputFile >> m_vectorTask2Priority[i];
			if (m_vectorTask2Priority[i] != -1) m_vectorPriority2Task[m_vectorTask2Priority[i]] = i;
		}
		cInputFile.close();
	}

	void PriorityAssignment::GenerateDkCPA(MultiProcessorSystem & rcSystem) {
		*this = PriorityAssignment(rcSystem);
		vector<int> vecSortedTasks(rcSystem.getNumTasks());
		for (int i = 0; i < rcSystem.getNumTasks(); i++) vecSortedTasks[i] = i;
		double m = rcSystem.getNumProcessors();
		double k = 0.5 * (m - 1 + sqrt(5 * pow(m, 2.0) - 6 * m + 1)) / m;
		auto score = [&](int iTaskIndex)->double {return rcSystem.getTask(iTaskIndex).D() - k * rcSystem.getTask(iTaskIndex).C();};
		sort(vecSortedTasks.begin(), vecSortedTasks.end(), [&](int lhs, int rhs)->bool { return score(lhs) < score(rhs); });
		for (int i = 0; i < rcSystem.getNumTasks(); i++) setPriority(vecSortedTasks[i], i);
	}

	// Assign priorities according to Deadline monotic (DMPO)
	void PriorityAssignment::GenerateDMPA(MultiProcessorSystem & rcSystem) {
		*this = PriorityAssignment(rcSystem);
		vector<int> vecSortedTasks(rcSystem.getNumTasks());
		for (int i = 0; i < rcSystem.getNumTasks(); i++) vecSortedTasks[i] = i;
		sort(vecSortedTasks.begin(), vecSortedTasks.end(), [&](int lhs, int rhs)->bool { return rcSystem.getTask(lhs).D() < rcSystem.getTask(rhs).D(); });
		for (int i = 0; i < rcSystem.getNumTasks(); i++) setPriority(vecSortedTasks[i], i);
	}

	// Added part: Assign priorities according to DCMPO (D-C)
	void PriorityAssignment::GenerateDCMPA(MultiProcessorSystem & rcSystem) {
		*this = PriorityAssignment(rcSystem);
		vector<int> vecSortedTasks(rcSystem.getNumTasks());
		for (int i = 0; i < rcSystem.getNumTasks(); i++) vecSortedTasks[i] = i;
		sort(vecSortedTasks.begin(), vecSortedTasks.end(), [&](int lhs, int rhs)->bool { return rcSystem.getTask(lhs).D()-rcSystem.getTask(lhs).C() < rcSystem.getTask(rhs).D()-rcSystem.getTask(rhs).C(); });
		for (int i = 0; i < rcSystem.getNumTasks(); i++) setPriority(vecSortedTasks[i], i);
	}

	// @return schedulability using HP-MUTER
	// 0 -> unschedulable
	// 1 -> schedulable after all priority levels are assigned
	// 2 -> early exit
	int PriorityAssignment::GenerateHP_MUTER(MultiProcessorSystem & rcSystem, float THRESH, bool DkC_if_not_MUTER, bool Lateness)
	{
		*this = PriorityAssignment(rcSystem);
		auto m = rcSystem.getNumProcessors();
		auto n = rcSystem.getNumTasks();

		// Compute utilization of the entire taskset
		float Total_Util = 0;
		for(int j = 0; j < n; ++j)
		{
			Task T = rcSystem.getTask(j);
			Total_Util += ((float)T.C()) / ((float)T.T());
		}
		
		#if DEBUG
		{
			cout << "Sort by DkC+Slack; " << THRESH << "; " << m << "; " << n << "; Util = " << Total_Util << "; Timeout = " << TIMEOUT << endl;
		}
		#endif 

		for(int k = n - 1; k >= 0; --k) // priority level, starting at the lowest
		{
			#if DEBUG_FINE
			{
				cout << "\n================================================\n";
				cout << "k = " << k << endl;
			}
			#endif

			// candidates with RT < D; ith element is the parameter of task i
			map<ValueType, float> sched_candidates;

			// ith element is the parameter of task i
			map<ValueType, float> all_candidates;

			// iterate through the tasks in rcSystem to compute RT of all candidates for priority k
			for(int i = 0; i < n; ++i)
			{
				#if DEBUG
				{
					cout << "\t\ti = " << i << endl;
				}
				#endif

				// skip if the task already has a priority assigned to it
				if(getPriority(i) != -1) 
				{
					continue; 
				}
				
				// if k = 0, there are no other unassigned tasks
				if(k > 0) 
				{
					ValueType RT_curr_task;

					PriorityAssignment cumulPA(rcSystem);
					MultiProcessorSystem subSystem(m);
					
					// mapping from task index in subSystem (hp) to task index in rcSystem (original taskset)
					vector<int> sub_to_rc;

					// All of the currently unassigned tasks have a higher priority compared to 'i'
					// compile the set hp(i) and calculate their utilization to determine which PA should be used
					float Util = 0;
					for(int j = 0; j < n; ++j)
					{
						// skip if it is the current task or if priority was already assigned
						if((j == i) || (getPriority(j) != -1)) continue;
						
						Task T = rcSystem.getTask(j);
						subSystem.AddTask(T);
						sub_to_rc.push_back(j);
						Util += ((float)T.C()) / ((float)T.T());
					}

					// Debug
					#if DEBUG_FINE
					{
						if(m != subSystem.getNumProcessors())
						{
							cout << "PROP_CHANGE\n";
						}
						if(subSystem.getNumTasks() != sub_to_rc.size())
						{
							cout << "subSystem " << subSystem.getNumTasks() << " and sub_to_rc " << sub_to_rc.size() << " mismatch\n";
						}
						if(k != subSystem.getNumTasks())
						{
							cout << "WRONG_k; k = " << k << "; num tasks in subSystem = " << subSystem.getNumTasks() << " \n";
						}
						if(k != sub_to_rc.size())
						{
							cout << "WRONG_vector k = " << k << "; num tasks in sub_to_rc = " << sub_to_rc.size() << " \n";
						}
					}
					#endif

					#if DEBUG
					{
						// Print for debug
						cout << "\t\tsub_to_rc\n";
						cout << "\t\tsubSystem_index\trcSystem_index\n";
						for(int subindex = 0; subindex < sub_to_rc.size(); ++subindex)
						{
							cout << "\t\t" << subindex << "\t\t\t" << sub_to_rc[subindex] << endl;
						}
						cout << endl;
					}
					#endif
					
					// Find suitable PA for hp(i)
					PriorityAssignment newPA;
					if(Util > (THRESH * (float)(m)))
					{
						#if DEBUG_FINE
						{
							cout << "\t\tUtil = " << Util << "; +++++++++DkC" << endl;
						}
						#endif
						newPA.GenerateDkCPA(subSystem);
					} 
					else
					{
						#if DEBUG_FINE
						{
							cout << "\t\tUtil = " << Util << "; &&&&&&&&&&&MUTER" << endl;
						}
						#endif
						MPSchedulingFR cFR(subSystem);
						cFR.setCorePerIter(5);
						auto MUTERPA = cFR.FocusRefinement(0, TIMEOUT);
						if(MUTERPA != 1) 
						{
							#if DEBUG
							{
								cout << "MUTER returned infeasible hp(i) at Util of hp(i) = " << Util << " for k = " << k 
							     << " and i = " << i << endl;
							}
							#endif

							if(DkC_if_not_MUTER)
							{
								PriorityAssignment newerPA;
								newerPA.GenerateDkCPA(subSystem);
								newPA = newerPA;
							}
							else
							{					
								// if MUTER returns non-feasible or timeout, then we skip the current task and move to the next one
								continue;
							}
						}
						else
						{
							newPA = cFR.getSolPA();
						}
					}
					
					#if DEBUG
					{
						// Print for debug
						cout << "\t\tNewPA (hp(i) after priority is assigned using MUTER or DkC; C, T, D should be noted for checking correctness instead of index)\n";
						cout << "\t\tsubindex\tpriority\tC\tD\tT\tR\trcIndex\t\tC\tD\tT\tR\n";
						for(int priority = 0; priority < sub_to_rc.size(); ++priority)
						{
							int index = newPA.getTask(priority);
							cout << "\t\t" << index << "\t\t" << priority;
							cout << "\t\t" << subSystem.getTask(index).C() << "\t" << subSystem.getTask(index).D() << "\t" << subSystem.getTask(index).T() << "\t" << subSystem.getTask(index).D();
							
							int rcIndex = sub_to_rc[index];
							cout << "\t" << rcIndex;
							cout << "\t\t" << rcSystem.getTask(rcIndex).C() << "\t" << rcSystem.getTask(rcIndex).D() << "\t" << rcSystem.getTask(rcIndex).T() << "\t" << rcSystem.getTask(rcIndex).D() << endl;
						}
						cout << endl;
					}
					#endif

					// add hp(i) to cumulPA
					for(int subIndex = 0; subIndex < sub_to_rc.size(); ++subIndex)
					{	
						// rcSystem index, subSystem priority
						cumulPA.setPriority(sub_to_rc[subIndex], newPA.getPriority(subIndex));

						#if DEBUG
						{
							cout << "\t\thp(i): Task " << sub_to_rc[subIndex] << " in rc assigned priority " 
								 << newPA.getPriority(subIndex) << endl;
						}
						// cout << endl;
						#endif
					}

					// add tasks whose priority has already been assigned to cumulPA
					for(int l = 0; l < n; ++l)
					{
						if(getPriority(l) == -1) continue; // skip unassigned tasks (which includes hp(i) and i)
						cumulPA.setPriority(l, getPriority(l));
						
						#if DEBUG
						{
							cout << "\t\talready assigned " << l << " to priority " << getPriority(l) << endl;
						}
						#endif
					}

					// assign k to the current task
					cumulPA.setPriority(i, k);

					#if DEBUG
					{
						// Print for debug
						cout << "\t\tcumulPA\n";
						cout << "\t\tindex\tpriority\tC\tD\tT\tR\n";
						for(int priority = 0; priority < n; ++priority)
						{
							int index = cumulPA.getTask(priority);
							cout << "\t\t" << index << "\t" << priority;
							cout << "\t\t" << rcSystem.getTask(index).C() << "\t" << rcSystem.getTask(index).D() << "\t" << rcSystem.getTask(index).T() << "\t" << rcSystem.getTask(index).D() << endl;
						}
						cout << endl;
					}
					#endif

					// Compute RT of task i at priority k
					RT_curr_task = ResponseTimeAnalysis_EPE().ComputeRTAGSYY2(rcSystem, cumulPA, k);

					#if DEBUG_FINE
					{
						cout << "\ti = " << i << "; D = " << rcSystem.getTask(i).D() << "; RT = " << RT_curr_task;
					}
					#endif

					// check if RT is <= D and if true, then add it to map
					if(RT_curr_task <= rcSystem.getTask(i).D())
					{
						// if estimation of hp(i) is correct, exit and return schedulable
						auto cEPE_WayOut = ResponseTimeAnalysis_EPE().ComputeRTAGSYY2(rcSystem, cumulPA);
						if(cEPE_WayOut.first)
						{
							#if DEBUG
							{
								cout << "EARLY_EXIT at k = " << k << "; i = " << i << "; result: " << cEPE_WayOut.first 
								 << "; Total_Util = " << Total_Util << "\n";
							}
							#endif
							return 2;
						}

						// if estimation of hp(i) is incorrect for the overall taskset, find DkC
						float C_sched = rcSystem.getTask(i).C();
						float D_sched = rcSystem.getTask(i).D();
						double m_sched = rcSystem.getNumProcessors();
						double k_sched = 0.5 * (m_sched - 1 + sqrt(5 * pow(m_sched, 2.0) - 6 * m_sched + 1)) / m_sched;
						sched_candidates[i] =  D_sched - (k_sched * C_sched);
						#if DEBUG_FINE
						{
							cout << "; DkC[i] = " << sched_candidates[i] << endl;
						}
						#endif
					}
					
					// if RT > D for all candidates, then slack is used to find the best candidate
					else
					{
						if(Lateness)
						{
							all_candidates[i] =  rcSystem.getTask(i).D() - RT_curr_task;
							#if DEBUG_FINE
							{
								cout << "; slack[i] = " << all_candidates[i] << endl;
							}
							#endif
						}
					}
				}
			
				if(k == 0)
				{
					setPriority(i, k);
					#if DEBUG_FINE
					{
						cout << "=============Task " << i << " assigned priority " << k << endl;
					}
					#endif
					break;
				}
			}

			if(!Lateness)
			{
				if(!all_candidates.empty())
				{
					#if DEBUG
					{
						cout << "lateness is false, but all_candidates is not empty\n";
					}
					#endif
				}
			}

			if(k == 0)
			{
				break;
			}

			// if both sched_candidates and all_candidates are empty, i.e., no candidates are schedulable by MUTER
			// taskset is not schedulable
			if((k > 0) && (sched_candidates.empty()) && (all_candidates.empty()))
			{
				#if DEBUG
				{
					cout << "Neither P1 or P2 could be computed; k = " << k << "; Total_Util = " << Total_Util << "\n";
				}
				#endif
				return 0;
			}
			
			// k > 0; but sched_candidates is non-empty implies that atleast one candidate has RT < D
			if((k > 0) && (!(sched_candidates.empty())))
			{
				#if DEBUG
				{
					cout << "P1 has been computed " << sched_candidates.size() << " times; k = " << k << "; Total_Util = " << Total_Util << "\n";
				}
				#endif
				#if DEBUG_FINE
				{
					cout << "Sched_Candidates:\n";
					cout << "Index\tParam\n";
					for(auto st_iter = sched_candidates.begin(); st_iter != sched_candidates.end(); ++st_iter)
					{
						cout << st_iter->first << "\t" << st_iter->second << endl;
					}
				}
				#endif
				
				// find the schedulable task with highest param 
				auto st_iter = sched_candidates.begin();
				int max_index = st_iter->first;
				float max_param = st_iter->second;
				
				#if DEBUG_FINE
				{
					cout << "initial max index: " << max_index << endl;
				}
				#endif
				for(st_iter = sched_candidates.begin(); st_iter != sched_candidates.end(); ++st_iter)
				{
					if(max_param < st_iter->second)
					{
						max_param = st_iter->second;
						max_index = st_iter->first;
					}
				}

				// assign priority k to the schedulable candidate with the largest param
				setPriority(max_index, k);
				#if DEBUG_FINE
				{
					cout << "=============Task " << max_index << " with parameter " << max_param << " assigned priority " << k << endl;
				}
				#endif
			}

			// k > 0; but sched_candidates is empty implies that none of the candidates have RT < D
			if((k > 0) && (sched_candidates.empty()) && (!(all_candidates.empty())))
			{
				#if DEBUG
				{
					cout << "P1 has been computed " << sched_candidates.size() 
					     << " times and P2 has been computed " << all_candidates.size() 
						 << " times; k = " << k << "; Total_Util = " << Total_Util << "\n";
				}
				#endif
				#if DEBUG_FINE
				{
					cout << "All_Candidates:\n";
					cout << "Index\tParam\n";
					for(auto st_iter = all_candidates.begin(); st_iter != all_candidates.end(); ++st_iter)
					{
						cout << st_iter->first << "\t" << st_iter->second << endl;
					}
				}
				#endif
				
				// find the task with max slack time
				auto st_iter = all_candidates.begin();
				int max_index = st_iter->first;
				float max_param = st_iter->second;
				
				#if DEBUG_FINE
				{
					cout << "initial max index: " << max_index << endl;
				}
				#endif
				for(st_iter = all_candidates.begin(); st_iter != all_candidates.end(); ++st_iter)
				{
					if(max_param < st_iter->second)
					{
						max_param = st_iter->second;
						max_index = st_iter->first;
					}
				}

				// assign priority k to the candidate with the largest slack time
				setPriority(max_index, k);
				#if DEBUG_FINE
				{
					cout << "=============Task " << max_index << " with parameter " << max_param << " assigned priority " << k << endl;
				}
				#endif
			}
		}

		// Test schedulability with EPE
		auto cEPETest = ResponseTimeAnalysis_EPE().ComputeRTAGSYY2(rcSystem, *this);
		#if DEBUG_FINE
		{
			cout << "==========================EPE result: " << cEPETest.first << endl;
		}
		#endif
		
		return cEPETest.first;
	}

	//RandomSystemGenerator
	RandomSystemGenerator::RandomSystemGenerator() {

	}

	vector<ValueType> RandomSystemGenerator::GeneratePeriod(int iTaskNum, ValueType lb, ValueType ub) {
		vector<ValueType> vectorPeriodTable(iTaskNum, 0);
		for (int i = 0; i < iTaskNum; i++) {
			vectorPeriodTable[i] = getRandomRealValue(lb, ub);
		}
		return vectorPeriodTable;
	}

	vector<ValueType> RandomSystemGenerator::GeneratePeriodLogUniform(int iTaskNum, ValueType lb, ValueType ub) {
		double dLogLB = log(lb), dLogUB = log(ub);
		vector<ValueType> vectorPeriodTable(iTaskNum, 0);
		for (int i = 0; i < iTaskNum; i++) {
			double dRand = getRandomRealValue(dLogLB, dLogUB);
			vectorPeriodTable[i] = round(exp(dRand));
		}
		return vectorPeriodTable;
	}

	vector<double> RandomSystemGenerator::UUnifastDiscard(int iTaskNum, double dTotalUtil) {
		bool bInvalid = true;
		vector<double> vectorUtil(iTaskNum, 0);
		while (bInvalid) {			
			double dSumU = dTotalUtil;
			for (int i = 0; i < iTaskNum - 1; i++){
				double dRand = getRandomRealValue(0.0, 1.0);
				double dNextSumU = dSumU * pow(dRand, (double)(1.0 / (iTaskNum - i - 1)));
				vectorUtil[i] = dSumU - dNextSumU;
				dSumU = dNextSumU;
				if ((bInvalid = (vectorUtil[i] > 1.0))) break;
			}
			vectorUtil.back() = (dSumU);
			bInvalid |= vectorUtil.back() > 1.0;
		}		
		return vectorUtil;
	}

	MultiProcessorSystem RandomSystemGenerator::GenerateSystem(int iTaskNum, int iProcessorNum, double dTotalUtil, ValueType cPeriodLB, ValueType cPeriodUB) {
		vector<double> vectorUtilTable = UUnifastDiscard(iTaskNum, dTotalUtil);		
		vector<ValueType> vectorPeriodTable = GeneratePeriodLogUniform(iTaskNum, cPeriodLB, cPeriodUB);		
		MultiProcessorSystem cSystem(iProcessorNum);
		for (int i = 0; i < iTaskNum; i++) {
			ValueType C = max(1.0, round((double)vectorPeriodTable[i] * (double)vectorUtilTable[i]));
			cSystem.AddTask(Task(C, vectorPeriodTable[i], vectorPeriodTable[i]));
		}
		return cSystem;
	}

	MultiProcessorSystem RandomSystemGenerator::GenerateSystemConstrainedDeadline(int iTaskNum, int iProcessorNum, double dTotalUtil, ValueType cPeriodLB, ValueType cPeriodUB) {
		vector<double> vectorUtilTable = UUnifastDiscard(iTaskNum, dTotalUtil);
		vector<ValueType> vectorPeriodTable = GeneratePeriodLogUniform(iTaskNum, cPeriodLB, cPeriodUB);
		MultiProcessorSystem cSystem(iProcessorNum);
		for (int i = 0; i < iTaskNum; i++) {
			ValueType C = max(1.0, round((double)vectorPeriodTable[i] * (double)vectorUtilTable[i]));
			ValueType D = getRandomIntValue(C, vectorPeriodTable[i]);
			cSystem.AddTask(Task(C, vectorPeriodTable[i], D));
		}
		return cSystem;
	}

	vector<double> RandomSystemGenerator::UUnifastMLforRT(int iTaskNum, double dTotalUtil) {
		bool bInvalid = true;
		vector<double> vectorUtil(iTaskNum, 0);

		// vectors of distributions
		vector<double> constants{0.1, 0.3, 0.5, 0.7, 0.9};

		// Random Number Generator Engine
		random_device rd;
		default_random_engine gen(rd());

		double constant = 0;
		int count = 0; // DEBUGGING
		
		while (bInvalid) {	
			count ++; //DEBUGGING
			// pick a random one
			ValueType dist = getRandomIntValue(0, 1); // 0 -> binomial, 1 -> exponential
			constant = constants[getRandomIntValue(0, 4)]; // pick random constant to use

			if(constant == 0.9) cout << "\tclosed interval\n"; //DEBUGGING
			
			// Distributions initialization
			//uniform_real_distribution<> uReal(0, 1);
			binomial_distribution<> binomial(iTaskNum, constant);
			exponential_distribution<> exponential(1 / constant);
			
			//double dSumU = dTotalUtil;
			double sumU = 0;
			for (int i = 0; i < iTaskNum; i++){
				
				double dRand = 0;
								
				if(dist == 0) // binomial
				{
					dRand = (double)(binomial(gen)) / (double)(iTaskNum);
				}
				else // exponential
				{
					dRand = exponential(gen);									
				}

				sumU += dRand; 
				//double dNextSumU = dSumU * pow(dRand, (double)(1.0 / (iTaskNum - i - 1)));
				vectorUtil[i] = dRand; //dSumU - dNextSumU;
				//dSumU = dNextSumU;
				if ((bInvalid = (sumU > dTotalUtil))) break;
			}
			//vectorUtil.back() = (dSumU);
			//bInvalid |= vectorUtil.back() > 1.0;
		}
		cout << "\tNo. of runs: " << count << endl; //DEBUGGING
		return vectorUtil;
	}

	MultiProcessorSystem RandomSystemGenerator::GenerateSystemMLforRT(int iTaskNum, int iProcessorNum, double dTotalUtil, ValueType cPeriodLB, ValueType cPeriodUB) {
		bool redo = false;
		
		// restart from here
		do
		{
			cout << "========Restart\n";
			vector<double> vectorUtilTable = UUnifastMLforRT(iTaskNum, dTotalUtil);	
			cout << "Utils generated\n";	
			vector<ValueType> vectorPeriodTable = GeneratePeriodLogUniform(iTaskNum, cPeriodLB, cPeriodUB);		
			MultiProcessorSystem cSystem(iProcessorNum);
			for (int i = 0; i < iTaskNum; i++) {
				ValueType C = max(1.0, round((double)vectorPeriodTable[i] * (double)vectorUtilTable[i]));
				cSystem.AddTask(Task(C, vectorPeriodTable[i], vectorPeriodTable[i]));
			}

			// RTA-LC = GSYY
			// If sched by any one of the following tests, regenerate the taskset
			
			// Test using RTA-LC + DkC
			PriorityAssignment RTALC_DkC; RTALC_DkC.GenerateDkCPA(cSystem);
			auto RTALC_DkC_Result = ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(cSystem, RTALC_DkC);
			if(redo = RTALC_DkC_Result.first)
			{ 
				cout << "Task set schedulable by RTA-LC + DkC done\n";
				continue;
			}
			cout << "Task set not schedulable by RTA-LC + DkC\n";
			
			// Test using RTA-LC + DCMPO
			PriorityAssignment RTALC_DCMPA; RTALC_DCMPA.GenerateDCMPA(cSystem);
			auto RTALC_DCMPA_Result = ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(cSystem, RTALC_DCMPA);
			if(redo = RTALC_DCMPA_Result.first)
			{ 
				cout << "Task set schedulable by RTA-LC + DCMPA\n";
				continue;
			}
			cout << "Task set not schedulable by RTA-LC + DCMPA\n";

			// Test using RTA-LC + DMPO
			PriorityAssignment RTALC_DMPA; RTALC_DMPA.GenerateDMPA(cSystem);
			auto RTALC_DMPA_Result = ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(cSystem, RTALC_DMPA);
			if(redo = RTALC_DMPA_Result.first)
			{ 
				cout << "Task set schedulable by RTA-LC + DMPA\n";
				continue;
			}
			cout << "Task set not schedulable by RTA-LC + DMPA\n";

			// Test sched using DA-LC + OPA
			PriorityAssignment DALC_OPA(cSystem);
			auto DALC_OPA_Result = DeadlineAnalysis().Audsley(cSystem, DALC_OPA, DA_LC);
			if(redo = DALC_OPA_Result)
			{
				cout << "Task set schedulable by DA-LC + OPA\n";
				continue;
			}
			cout << "Task set not schedulable by DA-LC + OPA\n";

			// Test sched using C-RTA + OPA
			PriorityAssignment CRTA_OPA(cSystem);
			auto CRTA_OPA_Result = DeadlineAnalysis().Audsley(cSystem, CRTA_OPA, C_RTA);
			if(redo = !CRTA_OPA_Result)
			{
				cout << "Task not schedulable by C-RTA + OPA\n";
				continue;
			}			
			cout << "Task schedulable by C-RTA + OPA\n";

			return cSystem;
		} while (redo);
	}
}