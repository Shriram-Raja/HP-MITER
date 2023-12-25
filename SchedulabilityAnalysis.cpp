#include "stdafx.h"
#include "SchedulabilityAnalysis.h"
#include "Util.h"
#include <algorithm>
#include <assert.h>
#include <set>
namespace MPScheduling {
	ResponseTimeAnalysis_Guan::ResponseTimeAnalysis_Guan() {
	
	}

	inline ValueType ResponseTimeAnalysis_Guan::CarryinInterference(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, ValueType dResponseTime) {
		auto C = rcSystem.getTask(iTaskIndex).C();
		auto T = rcSystem.getTask(iTaskIndex).T();
		return CLIP(((max((ValueType)0, t - C) % T) - (T - dResponseTime)), 0, C - 1);
	}

	inline ValueType ResponseTimeAnalysis_Guan::WCI(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, ValueType dResponseTime) {
		auto C = rcSystem.getTask(iTaskIndex).C();
		auto T = rcSystem.getTask(iTaskIndex).T();
		auto alpha = CarryinInterference(rcSystem, iTaskIndex, t, dResponseTime);
		return (MY_MAX(0, t - C) / T) * C + C + alpha;
	}

	inline ValueType ResponseTimeAnalysis_Guan::WNC(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t) {
		auto C = rcSystem.getTask(iTaskIndex).C();
		auto T = rcSystem.getTask(iTaskIndex).T();
		return (t / T)*C + MY_MIN(t % T, C);
	}
	
	ValueType ResponseTimeAnalysis_Guan::RHS(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA) {
		int iNumTasks = rcSystem.getNumTasks();
		int iThisPriority = rcPA.getPriority(iTaskIndex);
		int iNumProcessors = rcSystem.getNumProcessors();
		auto C = rcSystem.getTask(iTaskIndex).C();
		multiset<ValueType> heap;
		ValueType cWorkLoad = 0;		

		for (int i = 0; i < iNumTasks; i++) {
			if (rcPA.getPriority(i) >= iThisPriority) continue;
			auto nc = CLIP(WNC(rcSystem, i, t), 0, t - C + 1);
			auto ci = CLIP(WCI(rcSystem, i, t, rvectorResponseTime[i]), 0, t - C + 1);
			heap.insert(ci - nc);
			if (heap.size() > iNumProcessors - 1) heap.erase(heap.begin());
			cWorkLoad += nc;
		}
		for (auto ele : heap) cWorkLoad += ele;		
		return cWorkLoad / iNumProcessors + C;
	}

	ValueType ResponseTimeAnalysis_Guan::RT(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA) {
		ValueType t = rcSystem.getTask(iTaskIndex).C();
		ValueType t_RHS = RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA);
		while (t != t_RHS) {
			t = t_RHS;
			t_RHS = RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA);
		}
		return t;
	}

	double dRTTerminateOnBoundTime = 0; // not useful for now

	ValueType ResponseTimeAnalysis_Guan::RTTerminateOnBound(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA, ValueType cBound) {
		Timer cTimer; cTimer.Start();		
		ValueType t = rcSystem.getTask(iTaskIndex).C();
		ValueType t_RHS = RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA);
		while (t != t_RHS) {
			t = t_RHS;
			if (t > cBound) break;
			t_RHS = RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA);			
		}
		cTimer.Stop(); dRTTerminateOnBoundTime += cTimer.getElapsedTime_ms();
		return t;
	}

	ValueType ResponseTimeAnalysis_Guan::RT(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA, ValueType cLB, ValueType cUB) {
		Timer cTimer; cTimer.Start();
		ValueType t = max(cLB, rcSystem.getTask(iTaskIndex).C());
		ValueType t_RHS = RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA);
		while (t != t_RHS) {
			t = t_RHS;
			if (t > cUB) break;
			t_RHS = RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA);
		}
		cTimer.Stop(); dRTTerminateOnBoundTime += cTimer.getElapsedTime_ms();
		return t;
	}


	ValueType ResponseTimeAnalysis_Guan::RTNoCI(MultiProcessorSystem & rcSystem, int iTaskIndex, const PriorityAssignment & rcPA) {
		int iNumTasks = rcSystem.getNumTasks();
		int iThisPriority = rcPA.getPriority(iTaskIndex);
		int iNumProcessors = rcSystem.getNumProcessors();
		auto RHSNoCI = [&](ValueType t) {			
			auto C = rcSystem.getTask(iTaskIndex).C();			
			ValueType cWorkLoad = 0;
			for (int i = 0; i < iNumTasks; i++) {
				if (rcPA.getPriority(i) >= iThisPriority) continue;				
				cWorkLoad += WNC(rcSystem, i, t);
			}			
			return cWorkLoad / iNumProcessors + C;
		};
		ValueType t = rcSystem.getTask(iTaskIndex).C();		
		ValueType t_RHS = RHSNoCI(t);
		while (t != t_RHS) {
			t = t_RHS;
			t_RHS = RHSNoCI(t);
		}
		return t;
	}

	bool ResponseTimeAnalysis_Guan::ComputeAllRTGivenPA(MultiProcessorSystem & rcSystem, const PriorityAssignment & rcPA, vector<ValueType> & rvectorResponseTime){		
		assert(rvectorResponseTime.size() == rcSystem.getNumTasks());
		int iNumTask = rcSystem.getNumTasks();
		bool status = true;
		for (int iPriority = 0; iPriority < iNumTask; iPriority++) {
			int iTask = rcPA.getTask(iPriority);
			assert(iTask != -1);
			auto rt = RTTerminateOnBound(rcSystem, iTask, rvectorResponseTime, rcPA, rcSystem.getTask(iTask).D());
			rvectorResponseTime[iTask] = rt;
			status &= rt <= rcSystem.getTask(iTask).D();
		}
		return status;
	}

	pair<bool, vector<ValueType>> ResponseTimeAnalysis_Guan::ComputeAllRTGivenPA(MultiProcessorSystem & rcSystem, const PriorityAssignment & rcPA) {
		vector<ValueType> vectorRT(rcSystem.getNumTasks());
		bool bStatus = ComputeAllRTGivenPA(rcSystem, rcPA, vectorRT);
		return { bStatus, vectorRT };
	}

	// Implementation of C-RTA condition in Davis & Burns Tech Report 2010
	// C-RTA replaces the dResponseTime (RT of higher priority task i) with Ci thus giving the theoretical upper bound of RTA-LC (GSYY)
	// NOT A SCHED TEST. If a task set passes this, it could be schedulable by RTA-LC with some priority ordering, 
	// if not, it cannot be scheduled by RTA-LC in any priority ordering
	bool ResponseTimeAnalysis_Guan::CRTA(MultiProcessorSystem & rcSystem, int kTaskIndex, PriorityAssignment & rcPA){
		int iNumTask = rcSystem.getNumTasks();
		ValueType Ck = rcSystem.getTask(kTaskIndex).C(), Dk = rcSystem.getTask(kTaskIndex).D();
		if (Ck > Dk || Ck > rcSystem.getTask(kTaskIndex).T()) return false;

		// All hp(i) have RTi = Ci. Hence vector of Cs is created
		vector<ValueType> vectorRT(iNumTask);
		for (int i = 0; i < iNumTask; ++i)
		{
			vectorRT[i] = rcSystem.getTask(i).C();
		}

		auto result = RTTerminateOnBound(rcSystem, kTaskIndex, vectorRT, rcPA, Dk);
		return (result <= Dk);
	}


	// ****************************************************** Aug.5 Patch one: Add ZLL class for experiments ***********************************************

	ResponseTimeAnalysis_ZLL::ResponseTimeAnalysis_ZLL() {
	}

	// --------------------------------------------- Part 1 GSYY related computation ------------------------------------------------------
	inline ValueType ResponseTimeAnalysis_ZLL::CarryinInterference(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, ValueType dResponseTime) {
		auto C = rcSystem.getTask(iTaskIndex).C();
		auto T = rcSystem.getTask(iTaskIndex).T();
		return CLIP(((max((ValueType)0, t - C) % T) - (T - dResponseTime)), 0, C - 1);
	}

	inline ValueType ResponseTimeAnalysis_ZLL::WCI(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, ValueType dResponseTime) {
		auto C = rcSystem.getTask(iTaskIndex).C();
		auto T = rcSystem.getTask(iTaskIndex).T();
		auto alpha = CarryinInterference(rcSystem, iTaskIndex, t, dResponseTime);
		return (MY_MAX(0, t - C) / T) * C + C + alpha;
	}

	inline ValueType ResponseTimeAnalysis_ZLL::WNC(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t) {
		auto C = rcSystem.getTask(iTaskIndex).C();
		auto T = rcSystem.getTask(iTaskIndex).T();
		return (t / T)*C + MY_MIN(t % T, C);
	}
	
	ValueType ResponseTimeAnalysis_ZLL::RHS(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA) {
		int iNumTasks = rcSystem.getNumTasks();
		int iThisPriority = rcPA.getPriority(iTaskIndex);
		int iNumProcessors = rcSystem.getNumProcessors();
		auto C = rcSystem.getTask(iTaskIndex).C();
		multiset<ValueType> heap;
		ValueType cWorkLoad = 0;		

		for (int i = 0; i < iNumTasks; i++) {
			if (rcPA.getPriority(i) >= iThisPriority) continue;
			auto nc = CLIP(WNC(rcSystem, i, t), 0, t - C + 1);
			auto ci = CLIP(WCI(rcSystem, i, t, rvectorResponseTime[i]), 0, t - C + 1);
			heap.insert(ci - nc);
			if (heap.size() > iNumProcessors - 1) heap.erase(heap.begin());
			cWorkLoad += nc;
		}
		for (auto ele : heap) cWorkLoad += ele;		
		return cWorkLoad / iNumProcessors + C;
	}

	ValueType ResponseTimeAnalysis_ZLL::RT(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA) {
		ValueType t = rcSystem.getTask(iTaskIndex).C();
		ValueType t_RHS = RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA);
		while (t != t_RHS) {
			t = t_RHS;
			t_RHS = RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA);
		}
		return t;
	}


	ValueType ResponseTimeAnalysis_ZLL::RT(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA, ValueType cLB, ValueType cUB) {
		Timer cTimer; cTimer.Start();
		ValueType t = max(cLB, rcSystem.getTask(iTaskIndex).C());
		ValueType t_RHS = RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA);
		while (t != t_RHS) {
			t = t_RHS;
			if (t > cUB) break;
			t_RHS = RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA);
		}
		cTimer.Stop(); dRTTerminateOnBoundTime += cTimer.getElapsedTime_ms();
		return t;
	}


	ValueType ResponseTimeAnalysis_ZLL::RTNoCI(MultiProcessorSystem & rcSystem, int iTaskIndex, const PriorityAssignment & rcPA) {
		int iNumTasks = rcSystem.getNumTasks();
		int iThisPriority = rcPA.getPriority(iTaskIndex);
		int iNumProcessors = rcSystem.getNumProcessors();
		auto RHSNoCI = [&](ValueType t) {			
			auto C = rcSystem.getTask(iTaskIndex).C();			
			ValueType cWorkLoad = 0;
			for (int i = 0; i < iNumTasks; i++) {
				if (rcPA.getPriority(i) >= iThisPriority) continue;				
				cWorkLoad += WNC(rcSystem, i, t);
			}			
			return cWorkLoad / iNumProcessors + C;
		};
		ValueType t = rcSystem.getTask(iTaskIndex).C();		
		ValueType t_RHS = RHSNoCI(t);
		while (t != t_RHS) {
			t = t_RHS;
			t_RHS = RHSNoCI(t);
		}
		return t;
	}

	ValueType ResponseTimeAnalysis_ZLL::ComputeOmega(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA) {
		int iNumTasks = rcSystem.getNumTasks();
		int iThisPriority = rcPA.getPriority(iTaskIndex);
		int iNumProcessors = rcSystem.getNumProcessors();
		auto C = rcSystem.getTask(iTaskIndex).C();
		multiset<ValueType> heap;
		ValueType cWorkLoad = 0;		

		for (int i = 0; i < iNumTasks; i++) {
			if (rcPA.getPriority(i) >= iThisPriority) continue;
			auto nc = CLIP(WNC(rcSystem, i, t), 0, t - C + 1);
			auto ci = CLIP(WCI(rcSystem, i, t, rvectorResponseTime[i]), 0, t - C + 1);
			heap.insert(ci - nc);
			if (heap.size() > iNumProcessors - 1) heap.erase(heap.begin());
			cWorkLoad += nc;
		}
		for (auto ele : heap) cWorkLoad += ele;		
		return cWorkLoad;
	}

	// --------------------------------------------- Part 2 ZLL related computation ------------------------------------------------------

	ValueType ResponseTimeAnalysis_ZLL::ComputeMIN(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType x, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA)
	{   

		int iNumProcessors = rcSystem.getNumProcessors();
		auto Ci = rcSystem.getTask(iTaskIndex).C();

		int MIN_i_x = 0;
		int min_exec_time = 0;

		for(int y = 1; y <= x; y++)
		{
			min_exec_time =  std::min( std::max(y - ComputeOmega(rcSystem, iTaskIndex, y, rvectorResponseTime, rcPA) / iNumProcessors, 0) , Ci);
			if(min_exec_time > MIN_i_x)
			{
				MIN_i_x = min_exec_time;
			}
		}
		return MIN_i_x;
	}
	
	ValueType ResponseTimeAnalysis_ZLL::ComputeMAX(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType x, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA)
	{
		auto Ci = rcSystem.getTask(iTaskIndex).C();
		auto Ti = rcSystem.getTask(iTaskIndex).T();
		ValueType max_exec = 0;
		max_exec = Ci - ComputeMIN(rcSystem, iTaskIndex, Ti - x, rvectorResponseTime, rcPA);
		return max_exec;
	}

	ValueType ResponseTimeAnalysis_ZLL::ComputeBeta(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType x, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA)
	{
		auto Ci = rcSystem.getTask(iTaskIndex).C();
		auto Ti = rcSystem.getTask(iTaskIndex).T();
		ValueType gamma = std::max(x-Ci, 0) % Ti;
		ValueType calculated_MAX = ComputeMAX(rcSystem, iTaskIndex, gamma, rvectorResponseTime, rcPA);
		ValueType beta  = std::min( std::max(calculated_MAX,0), Ci - 1);
		return beta;
	}

	ValueType ResponseTimeAnalysis_ZLL::WCI_ZLL(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType x, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA){
		auto Ci = rcSystem.getTask(iTaskIndex).C();
		auto Ti = rcSystem.getTask(iTaskIndex).T();
		auto beta = ComputeBeta(rcSystem, iTaskIndex, x, rvectorResponseTime, rcPA); 
		return (std::max(0, x - Ci) / Ti) * Ci + Ci + beta;
	}


	ValueType ResponseTimeAnalysis_ZLL::RHS_ZLL(MultiProcessorSystem & rcSystem, int kTaskIndex, ValueType x, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA) {
		int iNumTasks = rcSystem.getNumTasks();
		int iThisPriority = rcPA.getPriority(kTaskIndex);
		int iNumProcessors = rcSystem.getNumProcessors();
		auto Ck = rcSystem.getTask(kTaskIndex).C();
		multiset<ValueType> heap;
		ValueType cWorkLoad = 0;		

		for (int i = 0; i < iNumTasks; i++) {
			if (rcPA.getPriority(i) >= iThisPriority) continue;
			auto nc = std::min( std::max(WNC(rcSystem, i, x), 0), x - Ck + 1 );
			auto ci = std::min( std::max( WCI_ZLL(rcSystem, i, x, rvectorResponseTime, rcPA), 0), x - Ck + 1);
			heap.insert(ci - nc);
			if (heap.size() > iNumProcessors - 1) heap.erase(heap.begin());
			cWorkLoad += nc;
		}
		for (auto ele : heap) cWorkLoad += ele;		
		return cWorkLoad / iNumProcessors + Ck;
	}

	ValueType ResponseTimeAnalysis_ZLL::RTTerminateOnBound(MultiProcessorSystem & rcSystem, int kTaskIndex, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA, ValueType cBound) {	
		ValueType t = rcSystem.getTask(kTaskIndex).C();
		ValueType t_RHS = RHS_ZLL(rcSystem, kTaskIndex, t, rvectorResponseTime, rcPA);
		while (t != t_RHS) {
			t = t_RHS;
			if (t > cBound) break;
			t_RHS = RHS_ZLL(rcSystem, kTaskIndex, t, rvectorResponseTime, rcPA);			
		}
		return t;
	}


	bool ResponseTimeAnalysis_ZLL::ComputeAllRTGivenPA(MultiProcessorSystem & rcSystem, const PriorityAssignment & rcPA, vector<ValueType> & rvectorResponseTime){		
		assert(rvectorResponseTime.size() == rcSystem.getNumTasks());
		int iNumTask = rcSystem.getNumTasks();
		bool status = true;

		// PHASE 1: Perform GSYY, if schedulable, return success
	    auto cDkCPAResult = ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(rcSystem, rcPA);
		if(cDkCPAResult.first != 0)
		{
			return true;
		}

		// ---------------- To calculate MAXi(x), it is based on rvectorResponse[iTask], if resp time doesnt change, MAXi wont change --
		// Therefore we consider response time instead of MAXi
		vector<ValueType> prev_responseTime = rvectorResponseTime;
		int maxIter = 10;
		int countIter = 0;
		int check_flag = 0;
		// -----------------------------------------------------------------

		do
	 	{ 
			for (int iPriority = 0; iPriority < iNumTask; iPriority++) {
				int iTask = rcPA.getTask(iPriority);
				assert(iTask != -1);
				auto rt = RTTerminateOnBound(rcSystem, iTask, rvectorResponseTime, rcPA, rcSystem.getTask(iTask).D());
				rvectorResponseTime[iTask] = rt;
				status &= rt <= rcSystem.getTask(iTask).D();
			}
			if(status==true) break;

			// ---------------- To calculate MAXi(x), it is based on rvectorResponse[iTask], if resp time doesnt change, MAXi wont change --
			vector<ValueType> curr_responseTime = rvectorResponseTime;
			if(curr_responseTime == prev_responseTime) check_flag=1;
			countIter += 1;
			if(countIter > maxIter)
			{
				status = false;
				break;
			}

	 	}while(check_flag == 0);

		return status;
	}

	pair<bool, vector<ValueType>> ResponseTimeAnalysis_ZLL::ComputeAllRTGivenPA(MultiProcessorSystem & rcSystem, const PriorityAssignment & rcPA) {
		vector<ValueType> vectorRT(rcSystem.getNumTasks());
		bool bStatus = ComputeAllRTGivenPA(rcSystem, rcPA, vectorRT);
		return { bStatus, vectorRT };
	}


	// *********************** Aug.5 Patch one: Add ZLL class for experiments ******************************

	// ---------------------------------- Part 1: GSYY base ------------------------------------------------
	// Constructor
	ResponseTimeAnalysis_EPE::ResponseTimeAnalysis_EPE() {

	}

	// Calculate alpha {below Eq. 2}
	inline ValueType ResponseTimeAnalysis_EPE::GSYY_CarryinInterference(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, ValueType dResponseTime){
		auto C = rcSystem.getTask(iTaskIndex).C();
		auto T = rcSystem.getTask(iTaskIndex).T();
		//cout << "\nAlpha: " << CLIP(((max((ValueType)0, t - C) % T) - (T - dResponseTime)), 0, C - 1);
		return CLIP(((max((ValueType)0, t - C) % T) - (T - dResponseTime)), 0, C - 1);
	}

	// Calculate Wk^CI {Eq. 1}
	inline ValueType ResponseTimeAnalysis_EPE::GSYY_WCI(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, ValueType dResponseTime) {
		auto C = rcSystem.getTask(iTaskIndex).C();
		auto T = rcSystem.getTask(iTaskIndex).T();
		auto alpha = GSYY_CarryinInterference(rcSystem, iTaskIndex, t, dResponseTime);
		//cout << "\nEq1; i = " << iTaskIndex << " Ci = " << C << " Ti = " << T << " Wci: " << (MY_MAX(0, t - C) / T) * C + C + alpha;
		return (MY_MAX(0, t - C) / T) * C + C + alpha;
	}

	// Calculate Wk^NC {Eq. 2}
	inline ValueType ResponseTimeAnalysis_EPE::GSYY_WNC(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t) {
		auto C = rcSystem.getTask(iTaskIndex).C();
		auto T = rcSystem.getTask(iTaskIndex).T();
		//cout << "\nEq2; i = " << iTaskIndex << " Ci = " << C << " Ti = " << T << " Wnc: " << (t / T)*C + MY_MIN(t % T, C);
		return (t / T)*C + MY_MIN(t % T, C);
	}
	
	// Calculate RHS of {Eq. 6}
	ValueType ResponseTimeAnalysis_EPE::GSYY_RHS(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, vector<vector<ValueType>> & rvectorResponseTime, const PriorityAssignment & rcPA) {
		int iNumTasks = rcSystem.getNumTasks();
		int iThisPriority = rcPA.getPriority(iTaskIndex);
		int iNumProcessors = rcSystem.getNumProcessors();
		auto C = rcSystem.getTask(iTaskIndex).C();
		multiset<ValueType> heap; // Store the interference {Eq.3 & 4} of higher priority tasks ascendingly, allow duplicate values
		ValueType cWorkLoad = 0;		

		for (int i = 0; i < iNumTasks; i++) {
			if (rcPA.getPriority(i) >= iThisPriority) continue; // The larger the number, the smaller the priority
			auto nc = CLIP(GSYY_WNC(rcSystem, i, t), 0, t - C + 1); 						// {Eq. 4}
			if(rvectorResponseTime[i].back() == rvectorResponseTime[i][-1]) cout << "-1 is correct";
			auto ci = CLIP(GSYY_WCI(rcSystem, i, t, rvectorResponseTime[i].back()/*[-1]*/), 0, t - C + 1); // {Eq. 3}
			//cout << "\nEq3; i = " << iTaskIndex << " Ici: " << ci;
			//cout << "\nEq4; i = " << iTaskIndex << " Ici: " << nc;
			heap.insert(ci - nc);
			if (heap.size() > iNumProcessors - 1){
				heap.erase(heap.begin());   	// {Eq. 10} The number of interference tasks is upper bounded by m - 1
			}
			cWorkLoad += nc; 					// Add all I_k^NC
		}
		for (auto ele : heap) cWorkLoad += ele;	// Add m - 1 largest I_k^CI
		//cout << "\nEq6 RHS: " << cWorkLoad / iNumProcessors + C;
		return cWorkLoad / iNumProcessors + C;  // {Eq. 6 RHS}
	}

	// Calculate minimum solution of Inequality {Eq. 6}
	ValueType ResponseTimeAnalysis_EPE::GSYY_RT(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<vector<ValueType>> & rvectorResponseTime, const PriorityAssignment & rcPA) {
		ValueType t = rcSystem.getTask(iTaskIndex).C();
		ValueType t_RHS = GSYY_RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA);
		while (t != t_RHS) { // Iterative fixed point search Eq. 6
			t = t_RHS;
			t_RHS = GSYY_RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA);
		}
		return t;
	}

	// Terminate Eq. 6 when response time exceeds the deadline (i.e.,cBound) of the task
	// double dRTTerminateOnBoundTime = 0; Already defined in L62
	ValueType ResponseTimeAnalysis_EPE::GSYY_RTTerminateOnBound(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<vector<ValueType>> & rvectorResponseTime, const PriorityAssignment & rcPA, ValueType cBound) {
		Timer cTimer; 
		cTimer.Start();		
		ValueType t = rcSystem.getTask(iTaskIndex).C();
		ValueType t_RHS = GSYY_RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA);
		while (t != t_RHS) {
			t = t_RHS;
			if (t > cBound) break;
			t_RHS = GSYY_RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA);			
		}
		cTimer.Stop(); 
		dRTTerminateOnBoundTime += cTimer.getElapsedTime_ms();
		return t;
	}

	// May not be useful if we do not provide lowerbound cLB
	ValueType ResponseTimeAnalysis_EPE::GSYY_RT(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<vector<ValueType>> & rvectorResponseTime, const PriorityAssignment & rcPA, ValueType cLB, ValueType cUB) {
		Timer cTimer; 
		cTimer.Start();
		ValueType t = max(cLB, rcSystem.getTask(iTaskIndex).C());
		ValueType t_RHS = GSYY_RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA);
		while (t != t_RHS) {
			t = t_RHS;
			if (t > cUB) break;
			t_RHS = GSYY_RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA);
		}
		cTimer.Stop(); 
		dRTTerminateOnBoundTime += cTimer.getElapsedTime_ms();
		return t;
	}

	// Calculate response time when the interference does NOT contain any carry-in jobs, currently it is not used
	ValueType ResponseTimeAnalysis_EPE::GSYY_RTNoCI(MultiProcessorSystem & rcSystem, int iTaskIndex, const PriorityAssignment & rcPA) {
		int iNumTasks = rcSystem.getNumTasks();
		int iThisPriority = rcPA.getPriority(iTaskIndex);
		int iNumProcessors = rcSystem.getNumProcessors();
		auto RHSNoCI = [&](ValueType t) {			
			auto C = rcSystem.getTask(iTaskIndex).C();			
			ValueType cWorkLoad = 0;
			for (int i = 0; i < iNumTasks; i++) {
				if (rcPA.getPriority(i) >= iThisPriority) continue;				
				cWorkLoad += GSYY_WNC(rcSystem, i, t);
			}			
			return cWorkLoad / iNumProcessors + C;
		};
		ValueType t = rcSystem.getTask(iTaskIndex).C();		
		ValueType t_RHS = RHSNoCI(t);
		while (t != t_RHS) {
			t = t_RHS;
			t_RHS = RHSNoCI(t);
		}
		return t;
	}

	// Final GSYY schedulability test to determine whether response time rt <= D
	bool ResponseTimeAnalysis_EPE::GSYY_ComputeAllRTGivenPA(MultiProcessorSystem & rcSystem, const PriorityAssignment & rcPA, vector<vector<ValueType>> & rvectorResponseTime){		
		assert(rvectorResponseTime.size() == rcSystem.getNumTasks());
		int iNumTask = rcSystem.getNumTasks();
		bool status = true;
		for (int iPriority = 0; iPriority < iNumTask; iPriority++) {
			int iTask = rcPA.getTask(iPriority);
			assert(iTask != -1);
			auto rt = GSYY_RTTerminateOnBound(rcSystem, iTask, rvectorResponseTime, rcPA, rcSystem.getTask(iTask).D());
			if(rvectorResponseTime[iTask].back() == rvectorResponseTime[iTask][-1]) cout << "-1 is correct";
			rvectorResponseTime[iTask].back()/*[-1]*/ = rt;
			status &= rt <= rcSystem.getTask(iTask).D(); // logic bit operation: true & false = false
		}
		return status;
	}

	// GSYY schedulability interface API
	pair<bool, vector<vector<ValueType>>> ResponseTimeAnalysis_EPE::GSYY_ComputeAllRTGivenPA(MultiProcessorSystem & rcSystem, const PriorityAssignment & rcPA) {
		vector<vector<ValueType>> vectorRT(rcSystem.getNumTasks(), vector<ValueType>(1, 0)); // Initialize with all 0s
		bool bStatus = GSYY_ComputeAllRTGivenPA(rcSystem, rcPA, vectorRT);
		return { bStatus, vectorRT };
		/*for(int i = 0; i < vectorRT.size(); i++){
			for(int j = 0; j < vectorRT[i].size(); j++){
				cout << "\nGSYY: Task " << i << " " << vectorRT[i][j];
			}
		}*/
	}

	// ---------------------------------- Part 2a: Overloaded GSYY for EPE ------------------------------------------------

	// Calculate RHS of {Eq. 6}
	ValueType ResponseTimeAnalysis_EPE::GSYY_RHS(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, vector<vector<ValueType>> & rvectorResponseTime, const PriorityAssignment & rcPA, int a) {
		int iNumTasks = rcSystem.getNumTasks();
		int iThisPriority = rcPA.getPriority(iTaskIndex);
		int iNumProcessors = rcSystem.getNumProcessors();
		auto C = a + 1; // a is the C_k value for which we are calculating RT
		//cout << "\nC = " << C << " x = " << t;
		multiset<ValueType> heap; // Store the interference {Eq.3 & 4} of higher priority tasks ascendingly, allow duplicate values
		ValueType cWorkLoad = 0;		

		for (int i = 0; i < iNumTasks; i++) {
			if (rcPA.getPriority(i) >= iThisPriority) continue; // The larger the number, the smaller the priority
			auto wnc = GSYY_WNC(rcSystem, i, t);
			if(rvectorResponseTime[i].back() == rvectorResponseTime[i][-1]) cout << "-1 is correct";
			auto wci = GSYY_WCI(rcSystem, i, t, rvectorResponseTime[i].back()/*[-1]*/);
			auto nc = CLIP(wnc, 0, t - C + 1); 						// {Eq. 4}
			auto ci = CLIP(wci, 0, t - C + 1); // {Eq. 3}
			//cout << "\nEq3; i = " << iTaskIndex << " Ici: " << ci;
			//cout << "\nEq4; i = " << iTaskIndex << " Inc: " << nc;
			heap.insert(ci - nc);
			//cout << "\nValue inserted into heap = " << ci - nc;
			if (heap.size() > iNumProcessors - 1){
				heap.erase(heap.begin());   	// {Eq. 10} The number of interference tasks is upper bounded by m - 1
			}
			cWorkLoad += nc; 					// Add all I_k^NC
			//cout << "\nChanging Workload = " << cWorkLoad;
		}
		for (auto ele : heap)
		{
			cWorkLoad += ele;	// Add m - 1 largest I_k^CI
			//cout << "\nEle = " << ele << " cWorkLoad = " << cWorkLoad;
		}
		//cout << "\nEq6 Omega: " << cWorkLoad << " RHS: " << (cWorkLoad / iNumProcessors) + C;
		return (cWorkLoad / iNumProcessors) + C;  // {Eq. 6 RHS}
	}

	// Terminate Eq. 6 when response time exceeds the deadline (i.e.,cBound) of the task
	ValueType ResponseTimeAnalysis_EPE::GSYY_RTTerminateOnBound(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<vector<ValueType>> & rvectorResponseTime, const PriorityAssignment & rcPA, ValueType cBound, int a) {
		Timer cTimer; 
		cTimer.Start();		
		ValueType t = a + 1; // current C value
		//cout << __FUNCTION__ << "x = " << t;
		ValueType t_RHS = GSYY_RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA, a);
		while (t != t_RHS) {
			//cout << "\n&&&&&&&&&&&&&&&&&&&&&&&Start&&&&&&&&&&&&&&&&&&&&&&&&&&&  GSYY  t: " << t << " RHS: " << t_RHS;
			t = t_RHS;
			if (t > cBound) break;
			t_RHS = GSYY_RHS(rcSystem, iTaskIndex, t, rvectorResponseTime, rcPA, a);

		}
		cTimer.Stop(); 
		dRTTerminateOnBoundTime += cTimer.getElapsedTime_ms();
		//cout << "\n&&&&&&&&&&&&&&&&&&&&&&&End&&&&&&&&&&&&&&&&&&&&&&&&&&&  GSYY  t: " << t << " RHS: " << t_RHS;
		return t;
	}


	// ---------------------------------- Part 2b: EPE related calculation ------------------------------------------------

	pair<bool, vector<ValueType>> ResponseTimeAnalysis_EPE::ComputeRTAGSYY2(MultiProcessorSystem & rcSystem, const PriorityAssignment & rcPA)
	{
		int n = rcSystem.getNumTasks();
		int m = rcSystem.getNumProcessors();

		vector<vector<ValueType>> taskset; // vector holding task parameters in order of priority
		// cout << "\t\t\tPrinting inside EPE2 interface\n\t\t\tindex\tpriority\tC\tD\tT\tR\n";
		for(int priority = 0; priority < n; ++priority)
		{
			int index = rcPA.getTask(priority);
			// cout << "\t\t\t" << index << "\t" << priority;
			// cout << "\t\t" << rcSystem.getTask(index).C() << "\t" << rcSystem.getTask(index).D() << "\t" << rcSystem.getTask(index).T() << "\t" << rcSystem.getTask(index).D() << endl;
			// C, D, T, worst-case R (filled in the GSYY function)
			taskset.push_back({rcSystem.getTask(index).C(), rcSystem.getTask(index).D(), rcSystem.getTask(index).T(), 0});
		}

		vector<ValueType> schedulability(n);

		int res = GSYY2(n, m, taskset, schedulability);

		bool result = (res != 0)? true: false;

		//cout << "\t\t---------------------------> GSYY2: " << res << endl;

		vector<ValueType> vectorRT; // vector of response times in index order
		for(int index = 0; index < n; ++index)
		{
			int priority = rcPA.getPriority(index);
			vectorRT.push_back(taskset[priority][3]);
		}
		return { result, vectorRT };		
	}

	// GSYY2 (EPE) interface to compute response time of task at priority 'curr_pr'
	ValueType ResponseTimeAnalysis_EPE::ComputeRTAGSYY2(MultiProcessorSystem & rcSystem, const PriorityAssignment & rcPA, int curr_pr)
	{
		//cout << "Entering interface\n";
		int n = rcSystem.getNumTasks();
		int m = rcSystem.getNumProcessors();

		vector<vector<ValueType>> taskset; // vector holding task parameters in order of priority
		//cout << "\t\t\tPrinting inside EPE2 interface\n\t\t\tindex\tpriority\tC\tD\tT\tR\n";
		for(int priority = 0; priority < n; ++priority)
		{
			int index = rcPA.getTask(priority);
			//cout << "\t\t\t" << index << "\t" << priority;
			// cout << "\t\t" << rcSystem.getTask(index).C() << "\t" << rcSystem.getTask(index).D() << "\t" << rcSystem.getTask(index).T() << "\t" << rcSystem.getTask(index).D() << endl;
			// C, D, T, worst-case R (filled in the GSYY function)
			taskset.push_back({rcSystem.getTask(index).C(), rcSystem.getTask(index).D(), rcSystem.getTask(index).T(), 0});
		}

		vector<ValueType> schedulability(n);

		ValueType responseTime = GSYY2(n, m, taskset, schedulability, curr_pr);

		//cout << "At the interface, response time: " << responseTime << endl;

		// vector<ValueType> vectorRT; // vector of response times in index order
		// for(int index = 0; index < n; ++index)
		// {
		// 	int priority = rcPA.getPriority(index);
		// 	vectorRT.push_back(taskset[priority][3]);
		// }
		return responseTime;
	}


	// ----------------------------- End of EPE class -------------------------------

	// Audsley's OPA
	// Sched Test: 
	// 0 - DA
	// 1 - DA-LC
	// 2 - C-RTA
	bool DeadlineAnalysis::Audsley(MultiProcessorSystem & rcSystem, PriorityAssignment & rcPA, int sched_test = 0) {
		int iTaskNum = rcSystem.getNumTasks();
		int iPriority = iTaskNum - 1;
		for (; iPriority >= 0; iPriority--)
		{
			bool bAssigned = false;
			for (int i = 0; i < iTaskNum && bAssigned == false; i++){
				if ((rcPA.getPriorityByTask(i) != -1)) continue; // skip already assigned tasks				
				rcPA.setPriority(i, iPriority); // assign current priority, iPriority to unassigned task, i
				
				// Test schedulability of task i at iPriority
				switch (sched_test)
				{
					case DA:
						bAssigned = isSchedulable(rcSystem, i, rcPA);
						break;

					case DA_LC:
						bAssigned = isSchedulable_LC(rcSystem, i, rcPA);
						break;

					case C_RTA:
						cout << "C-RTA + OPA\n";
						bAssigned = ResponseTimeAnalysis_Guan().CRTA(rcSystem, i, rcPA);
						break;
				
					default:
						cout << "Invalid test specified for OPA\n";
						break;
				}

				// If i in unschedulable at iPriority, unset priority of i and move to the next unassigned task
				if (!bAssigned) rcPA.unset(i);
			}
			// If all unassigned tasks have been traversed an none of them are schedulable at iPriority, then exit
			if (!bAssigned) break;
		}		
		return iPriority == -1;// if it equals -1, it means all the tasks got assigned a priority = successful. Otherwise OPA failed
	}

	// This is not used in our code.
	ValueType DeadlineAnalysis::Interference(MultiProcessorSystem & rcSystem, int iTaskIndex, PriorityAssignment & rcPA) {
		ValueType Ck = rcSystem.getTask(iTaskIndex).C(), Dk = rcSystem.getTask(iTaskIndex).D();
		if (Ck > Dk) return Dk + 1;// this should be return false
		ValueType cInterference = 0;
		for (int i = 0; i < rcSystem.getNumTasks(); i++) {
			if (rcPA.getPriority(i) >= rcPA.getPriority(iTaskIndex)) continue;
			cInterference += IDi(rcSystem, i, Dk, Ck);
		}
		cInterference = cInterference / rcSystem.getNumProcessors() + Ck;
		return cInterference;
	}

	// DeadlineAnalysis
	bool DeadlineAnalysis::isSchedulable(MultiProcessorSystem & rcSystem, int iTaskIndex, PriorityAssignment & rcPA) {		
		ValueType Ck = rcSystem.getTask(iTaskIndex).C(), Dk = rcSystem.getTask(iTaskIndex).D();
		if (Ck > Dk || Ck > rcSystem.getTask(iTaskIndex).T()) return false;
		ValueType cInterference = 0;
		for (int i = 0; i < rcSystem.getNumTasks(); i++) {
			if (rcPA.getPriority(i) >= rcPA.getPriority(iTaskIndex)) continue;							
			cInterference += IDi(rcSystem, i, Dk, Ck);						
		}
		cInterference = cInterference / rcSystem.getNumProcessors() + Ck; // Equation 2
		return cInterference <= Dk;
	}

	// Below Eqn.2. Review Robert Davis, exactly the same as equations
	ValueType DeadlineAnalysis::IDi(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType Dk, ValueType Ck){
		return min(WDi(rcSystem, iTaskIndex, Dk), Dk - Ck + 1);
	}

	ValueType DeadlineAnalysis::WDi(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t) {
		ValueType Ci = rcSystem.getTask(iTaskIndex).C(), Di = rcSystem.getTask(iTaskIndex).D(), Ti = rcSystem.getTask(iTaskIndex).T();
		return NDi(rcSystem, iTaskIndex, t) * Ci + min(Ci, t + Di - Ci - NDi(rcSystem, iTaskIndex, t) * Ti);
	}

	ValueType DeadlineAnalysis::NDi(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t) {
		ValueType Ci = rcSystem.getTask(iTaskIndex).C(), Di = rcSystem.getTask(iTaskIndex).D(), Ti = rcSystem.getTask(iTaskIndex).T();
		return (t + Di - Ci) / Ti;
	}


	// DeadlineAnalysis + Limited Carry-in
	bool DeadlineAnalysis::isSchedulable_LC(MultiProcessorSystem & rcSystem, int iTaskIndex, PriorityAssignment & rcPA) {		
		ValueType Ck = rcSystem.getTask(iTaskIndex).C(), Dk = rcSystem.getTask(iTaskIndex).D();
		if (Ck > Dk || Ck > rcSystem.getTask(iTaskIndex).T()) return false;
		ValueType cInterference = 0;
		vector<ValueType> IDIFF_D;
		for (int i = 0; i < rcSystem.getNumTasks(); i++) {
			if (rcPA.getPriority(i) >= rcPA.getPriority(iTaskIndex)) continue;
			auto nci = INCi(rcSystem, i, Dk, Ck);
			auto di = IDi(rcSystem, i, Dk, Ck);
			cInterference += nci;	
			IDIFF_D.push_back(di - nci);
		}
		sort(IDIFF_D.begin(), IDIFF_D.end(), greater<ValueType>());
		int carry_in = min((int)(IDIFF_D.size()), rcSystem.getNumProcessors() - 1);
		for (int lc_iter = 0; lc_iter < carry_in; lc_iter++) {
			cInterference += IDIFF_D[lc_iter];
		}
		cInterference = cInterference / rcSystem.getNumProcessors() + Ck; // Equation 14
		return cInterference <= Dk;
	}

	// Eqn.9 in DA-LC paper (Improved priority assignment for gFP pre-emptive sched in multiprocessor RTS, Davis & Burns)
	ValueType DeadlineAnalysis::INCi(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, ValueType Ck){
		return min(WNCi(rcSystem, iTaskIndex, t), t - Ck + 1);
	}

	ValueType DeadlineAnalysis::WNCi(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t) {
		ValueType Ci = rcSystem.getTask(iTaskIndex).C(), Di = rcSystem.getTask(iTaskIndex).D(), Ti = rcSystem.getTask(iTaskIndex).T();
		return (t / Ti) * Ci + MY_MIN(t % Ti, Ci);
	}

}