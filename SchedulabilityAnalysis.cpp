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


	// ---------------------------------- Add EPE class for experiments ------------------------------------
	// EPE is based on existing schedulability test GSYY
	// Notation: 
	// - CLIP is defined as [[a]]_b^c {Page 4, Section 2.1}
	// - t is the response time x in equation

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


	// ---------------------------------- Part 2: EPE related calculation ------------------------------------------------

	int ResponseTimeAnalysis_EPE::GSYY2(int n, int m, vector<vector<ValueType>> &taskset, vector<ValueType> &schedulability)
{
	double beginning, end, sumi, diffi[200];
	int flagchange, result, i, j, k, t, l, o, x, y, z, y1, y2, x1, x2, x3, x4, limitednumber, maxaddedwnc, maxaddedwci;
	int interference;
	int ce, ce1, ce2, ce3, ce4, ce5, ce6, ce7, ce8, ce9;
	flagchange = 1;
	int diffw[1000], wnc[1000], wci[1000];
	int wcrt[200][401];
	int ici[200], inc[200];
	int oici[200], oinc[200];
	int addedwnc[200], addedwci[200], addedinc[200], addedici[200];
	int minworkload[401];
	int wcet;
	k = 0;
	////////////////////////////////////////////////////////////////////////Ϊ���񸳳�ʼWCRT
	while (k<m)
	{
		taskset[k][3] = taskset[k][0];
		k++;
	}
	while (k<n)
	{
		taskset[k][3] = taskset[k][1];
		k++;
	}
	/////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////��ʼ������ɵ�����
	k = 0;
	while (k<m)
	{
		i = 1;
		while (i <= taskset[k][0])
		{
			wcrt[k][i] = i;
			i++;
		}
		schedulability[k] = 1;
		k++;
	}
	while (k<n)
	{
		schedulability[k] = 0;
		k++;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////��ʼ��wcrt
	k = 0;
	while (k<n)
	{
		wcrt[k][0] = 0;
		k++;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////
	k = m;
	while (k<n)
	{
		//////////////////////////////////////////////////////////////////////////��ʼ��inc��oici
		i = 0;
		while (i<n)
		{
			oinc[i] = 0;
			i++;
		}
		i = 0;
		while (i<n)
		{
			oici[i] = 0;
			i++;
		}
		i = 0;
		while (i < 401)
		{
			minworkload[i] = 0;
			i++;
		}
		//////////////////////////////////////////////////////////////////////////////////////////
		wcet = 1;
/*		for(int iteration = 0; iteration < 3; iteration++)
		{
			cout << "\nTask " << iteration << "\n";
			for(int jteration = 0; jteration < 4; jteration++)
			{
				cout << "Index: " << jteration << " value: " << taskset[iteration][jteration] << "\n";
			}
		}*/
		/////////////////////////////////////////////////////////////////////����Ck=wcetʱk��WCRT
		while (wcet <= taskset[k][0])
		{
			t = wcrt[k][wcet - 1] + 1;
			while (t <= taskset[k][1])
			{
				x = t - wcrt[k][wcet - 1];              //���䳤��
				i = 0;
				while (i<k)////////////////////////////////////////////////////////����i��k�ĸ���
				{
					//////////////////////////////////////////////////////////����t�����ڵ�W^{NC}
					if (t%taskset[i][2]>taskset[i][0])
					{
						wnc[i] = int(t / taskset[i][2])*taskset[i][0] + taskset[i][0];
					}
					else
					{
						wnc[i] = int(t / taskset[i][2])*taskset[i][0] + int(t%taskset[i][2]);
					}
					/////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////////////////////////////����t�����ڵ�W^{CI}
					if (t <= taskset[i][0])
					{
						wci[i] = t;
					}
					else
					{
						if ((t - taskset[i][0]) % taskset[i][2] > taskset[i][2] - taskset[i][3])
						{
							wci[i] = (int((t - taskset[i][0]) / taskset[i][2]) + 1)*taskset[i][0] + (t - taskset[i][0]) % taskset[i][2] - taskset[i][2] + taskset[i][3];
						}
						else
						{
							wci[i] = (int((t - taskset[i][0]) / taskset[i][2]) + 1)*taskset[i][0];
						}
					}
					/////////////////////////////////////////////////////////////////////////////////
					/////////////////////////////////////////////////////����WNC��WCI���������
					y = 0;
					maxaddedwnc = 0;
					while (y < taskset[i][2]&&y<=wcrt[k][wcet-1])
					{
						/////////////yʱ[0��wcrt[k][wcet-1])������i���nc����
						x1=(wcrt[k][wcet - 1] - y) / taskset[i][2] * taskset[i][0] + min(y, taskset[i][0]); ///////////////����i��[0��wcrt[k][wcet-1])�ϵ����nc���� 
						y1 = 0;///////////////////////////////////////////////����k��[wcrt[k][wcet-1]-y,wcrt[k][wcet-1])����Сִ��ʱ�䣻���y��Ӧ�Ĺ���������ֵ����,����i��[wcrt[k][wcet-1]-y,wcrt[k][wcet-1])��������k����С����
						while (wcrt[k][y1]+1 <= x + y)
						{
							y1++;
							if (y1 == wcet)
							{
								break;
							}
						}
						y1 = y1 - 1;
						if (y1 > y)
						{
							y++;
							continue;
						}
						if (y1 > taskset[i][0])
							y1 = taskset[i][0];
						if (y<= taskset[i][0])////////////////////////////////////////////////���x1������i��[0��wcrt[k][wcet-1])�϶�����k��������
						{
							x1 = x1 - y1;
						}
						else
						{
							y2 = 0;
							while (wcrt[k][y2] <= wcrt[k][wcet - 1] - y + taskset[i][0])
							{
								y2++;
							}
							y2 = y2 - 1;   
							x1 = x1-max(y1 + y2 - (wcet - 1), 0);
						}
						if (x1 < minworkload[i])
						{
							y++;
							continue;///////////////////////////////////wcrt[k][wcet-1]-y������Ϊ�ͷ�ʱ��
						}
						else
						{
							if (minworkload[i] <= (wcrt[k][wcet - 1] - y) / taskset[i][2] * taskset[i][0])/////////////y����Ӧ��Job����ִ��
							{
								if (taskset[i][3] <= y)
								{
									x2 = max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max((x + y - taskset[i][2]), 0) % taskset[i][2], taskset[i][0]);////x2ΪW�������
								}
								else
								{
									x2 = min(min(taskset[i][3]-y,taskset[i][0]-y1),x)+ max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max((x + y - taskset[i][2]), 0) % taskset[i][2], taskset[i][0]);
								}
							}
							else////////////////////////////////////////////////////[0,wcrt[k][wcet-1])��i����ִ��minworkload[i]-(wcrt[k][wcet - 1] - y) / taskset[i][2] * taskset[i][0]
							{
								y2 = minworkload[i] - (wcrt[k][wcet - 1] - y) / taskset[i][2] * taskset[i][0]+y1;////////////////////���y��Ӧ�Ĺ���������ֵ���أ���ô�˹�����[wcrt[k][wcet-1]-y,wcrt[k][wcet-1])�ϵ���С����
								if (taskset[i][3] <= y)
								{
									x2 = max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max((x + y - taskset[i][2]), 0) % taskset[i][2], taskset[i][0]);////x2ΪW�������
								}
								else
								{
									x2 = min(min(taskset[i][3] - y, max(taskset[i][0] - y2,0)), x) + max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max((x + y - taskset[i][2]), 0) % taskset[i][2], taskset[i][0]);
								}							
							}
							if (maxaddedwnc < x2)
							{
								maxaddedwnc = x2;
							}
						}
						y++;
					}
					y = 0;
					maxaddedwci = 0;
					while (y < taskset[i][2])
					{
						/////////////yʱ[0��wcrt[k][wcet-1])������i���nc����
						if (y <= wcrt[k][wcet - 1])
						{
							x1 = (wcrt[k][wcet - 1] - y) / taskset[i][2] * taskset[i][0] + min(y, taskset[i][0]) + max(min((wcrt[k][wcet - 1] - y) % taskset[i][2] - taskset[i][2] + taskset[i][3], taskset[i][0]), 0);
							y1 = 0;/////////////////////////    [ R_k^{wcet-1}-y, R_k^{wcet-1} )�ϵ���С����ʱ��  
							while (wcrt[k][y1]+1 <= x + y)
							{
								y1++;
								if (y1 == wcet)
								{
									break;
								}
							}
							y1 = y1 - 1;
							if (y1 > y)
							{
								y++;
								continue;
							}
							if (y1 > taskset[i][0])
								y1 = taskset[i][0];
							if (y <= taskset[i][0])
							{
								x1 = x1 - y1;
							}
							else
							{
								y2 = 0;
								while (wcrt[k][y2] <= wcrt[k][wcet - 1] - y + taskset[i][0])
								{
									y2++;
								}
								y2 = y2 - 1;
								x1 = x1 - max(y1 + y2 - (wcet - 1), 0);
							}
							if (x1 < minworkload[i])
							{
								y++;
								continue;///////////////////////////////////wcrt[k][wcet-1]-y������Ϊ�ͷ�ʱ��
							}
							else
							{
								if (minworkload[i] <= (wcrt[k][wcet - 1] - y) / taskset[i][2] * taskset[i][0] + max(min((wcrt[k][wcet - 1] - y) % taskset[i][2] - taskset[i][2] + taskset[i][3], taskset[i][0]), 0))/////////////y����Ӧ��Job����ִ��
								{
									if (taskset[i][3] <= y)
									{
										x2 = max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max(x + y - taskset[i][2], 0) % taskset[i][2], taskset[i][0]);
									}
									else
									{
										x2 = min(min(taskset[i][0]-y1, taskset[i][3]-y),x) + max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max(x + y - taskset[i][2], 0) % taskset[i][2], taskset[i][0]);
									}
								}
								else////////////////////////////////////////////////////[0,wcrt[k][wcet-1])��i����ִ��minworkload[i]-(wcrt[k][wcet - 1] - y) / taskset[i][2] * taskset[i][0]
								{
									y2 = minworkload[i] - (wcrt[k][wcet - 1] - y) / taskset[i][2] * taskset[i][0] - max(min((wcrt[k][wcet - 1] - y) % taskset[i][2] - taskset[i][2] + taskset[i][3], taskset[i][0]), 0)+y1;
									if (taskset[i][3] <= y)
									{
										x2 = max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max(x + y - taskset[i][2], 0) % taskset[i][2], taskset[i][0]);
									}
									else
									{
										x2 = min(min(max(taskset[i][0] - y2,0), taskset[i][3] - y),x) + max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max(x + y - taskset[i][2], 0) % taskset[i][2], taskset[i][0]);
									}
								}
								if (maxaddedwci < x2)
								{
									maxaddedwci = x2;
								}
							}
						}
						else///////////////////////y>wcrt[k][wcet-1]
						{
							y1 = wcet-1;///////////////////////////////////////////////////////////////�������i��wcrt[k][wcet-1]��ִ�У�����i��[0,wcrt[k][wcet-1])�ϵ���С����Ϊy1
							if (y1>taskset[i][0])
								y1 = taskset[i][0];
							if (taskset[i][3] >= y)
							{
								x1 = min(taskset[i][0], wcrt[k][wcet-1]);////////////////[0,wcrt[k][wcet-1])�ϵ������
								if (wcrt[k][wcet - 1] <= taskset[i][0])
									x1 = x1 - y1;////////////////////////////////////////[0,wcrt[k][wcet-1])�ϵ�������
								else
								{
									y2 = 0;
									while (wcrt[k][y2] <= taskset[i][0])
									{
										y2++;
										if (y2 == wcet)
										{
											break;
										}
									}
									y2 = y2 - 1;
									x1 = x1 - y2;
								}
							}
							else
							{
								x1 = max(min(taskset[i][0], taskset[i][3] + wcrt[k][wcet - 1] - y), 0);
								y2 = 0;
								while (wcrt[k][y2] <= x1)
								{
									y2++;
									if (y2 == wcet)
									{
										break;
									}
								}
								y2 = y2 - 1;
								x1 = x1 - y2;
							}
							if (x1 < minworkload[i])
							{
								y++;
								continue;///////////////////////////////////wcrt[k][wcet-1]-y������Ϊ�ͷ�ʱ��
							}
							else
							{
								y2 = minworkload[i] + y1;////////     �������i��wcrt[k][wcet-1]֮��ִ�У� ��Ҫ����+��С����
								if (taskset[i][3] <= y)
									x2 = max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max((x + y - taskset[i][2]), 0) % taskset[i][2], taskset[i][0]);        ////x2ΪW�������
								else
									x2 = min(min(max(min(taskset[i][0], taskset[i][3] - y + wcrt[k][wcet - 1]) - y2, 0),taskset[i][3]-y), x)+max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max((x + y - taskset[i][2]), 0) % taskset[i][2], taskset[i][0]);
								if (maxaddedwci < x2)
								{
									maxaddedwci = x2;
								}
							}
						}
						y++;
					}
					/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////�������W=I1+W2
					if (wnc[i]>oinc[i] + maxaddedwnc)
						wnc[i] = oinc[i] + maxaddedwnc;
					if (wci[i] > oici[i] + maxaddedwci)
						wci[i] = oici[i] + maxaddedwci;
					if (wnc[i] < oinc[i] + maxaddedwnc)
						wnc[i] = wnc[i] ;
					if (wci[i] < oici[i] + maxaddedwci)
						wnc[i] = wnc[i] ;
					/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					////////////////////////////////////////////////////////////////����t�����ڵ�I
					ce = t - wcet + 1;
					if (wci[i] > t - wcet + 1)
					{
						ici[i] = t - wcet + 1;
				//		limitednumber = limitednumber + 1;
					}
					else
					{
						ici[i] = wci[i];
					}
					if (wnc[i] > t - wcet + 1)
					{
						inc[i] = t - wcet + 1;
				//		limitednumber = limitednumber + 1;
					}
					else
					{
						inc[i] = wnc[i];
					}
					//////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////////////////////////////��¼i��t�ϵ�������
					//////////////////////////////////////////////////////////////////////////////
					diffw[i] = ici[i] - inc[i];///////////////////////////////////////����I^{diff}
					i++;
				}////////////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////�����ܸ���
				i = 0;
				interference = 0;
				while (i<k)
				{
					interference = interference + inc[i];
					i++;
				}
				j = 0;
				while (j<m - 1)
				{
					l = 0;
					o = 0;
					while (l<k)
					{
						if (diffw[o]<diffw[l])
						{
							o = l;
						}
						l++;
					}
					interference = interference + diffw[o];
					diffw[o] = 0;
					j++;
				}
				////////////////////////////////////////////////////////////////////////////////
				interference = int(interference / m);
				if (interference + wcet <= t)
				{
					wcrt[k][wcet] = t;
					/////////////////////////////////////////////////////����ÿ��������t�ϵ���С����
					j = 0;
					while (j < k)
					{
						if (ici[j]>interference)
						{
							ici[j] = interference;
						}
						if (inc[j]>interference)
							inc[j] = interference;
						j++;
					}
					i = 0;
					while (i < k)
					{
						j = 0;
						sumi= 0;
						while (j<k)
						{
							if (j != i)
							{
								sumi = sumi+ inc[i];
							}							
							j++;
						}
						j = 0;
						while (j < k)
						{
							diffi[j] = max(ici[j] - inc[j],0);
							j++;
						}
						j = 0;
						while (j<m - 1)
						{
							l = 0;
							o = 0;
							while (l<k)
							{
								if (l != i)
								{
									if (diffw[o] < diffw[l])
									{
										o = l;
									}
								}
								l++;
							}
							sumi = sumi + diffi[o];
							diffi[o] = 0;
							j++;
						}
						minworkload[i] = MY_MAX(interference*m - sumi,0);
						i++;
					}
					/////////////////////////////////////////////////////////////////////////////////
					if (wcet == taskset[k][0])
					{
						schedulability[k] = 1;
						taskset[k][3] = t;
					}
					break;
				}
				else
				{
					t = interference + wcet;
				}
			}
			if (t>taskset[k][1])  ////////////////////Rk>Dk
			{
				wcrt[k][wcet] = taskset[k][1];
				taskset[k][3] = taskset[k][1];
				schedulability[k] = 0;
				cout << "unsched at k = " << k << " because at a = " << wcet << " the wcrt = " << t << ", but D = " << taskset[k][1] << endl;
				break;
			}
			i = 0;
			while (i < k)
			{
				if (inc[i] <= t - wcet)
					oinc[i] = inc[i];
				else
					oinc[i] = t - wcet;
				if (oici[i] <= t - wcet)
					oici[i] = ici[i];
				else
					oici[i] = t-wcet;
				i++;
			}
			//cout << "k = " << k << " t = " << t << " wcet = " << wcet << "\n";
			wcet++;
		}
		////////////////////////////////////////////////////////////////////////////////////////
		//cout << "Outside: k = " << k << " t = " << t << " wcet = " << wcet << "\n";
		if (schedulability[k] == 0)
		{
			break;
		}
		////////////////////////////////////////////////////////////////////////////////////////
		k++;
	}
	//////////////////////////////////////////////////////////////////////////////////�жϿɵ�����
	i = 0;
	result = 1;
	while (i<n)
	{
		if (schedulability[i] == 0)
		{
			result = 0;
			break;
		}
		i++;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////
	return result;
}

// GSYY2 overloaded to return only response time of task at priority 'curr_pr'
	int ResponseTimeAnalysis_EPE::GSYY2(int n, int m, vector<vector<ValueType>> &taskset, vector<ValueType> &schedulability, int curr_pr)
{
	// if curr_pr is less than m, wcrt = wcet
	//cout << "curr_pr = " << curr_pr << endl;
	if(curr_pr < m)
	{
		return taskset[curr_pr][0];
	}
	double beginning, end, sumi, diffi[200];
	int flagchange, result, i, j, k, t, l, o, x, y, z, y1, y2, x1, x2, x3, x4, limitednumber, maxaddedwnc, maxaddedwci;
	int interference;
	int ce, ce1, ce2, ce3, ce4, ce5, ce6, ce7, ce8, ce9;
	flagchange = 1;
	int diffw[1000], wnc[1000], wci[1000];
	int wcrt[200][401];
	int ici[200], inc[200];
	int oici[200], oinc[200];
	int addedwnc[200], addedwci[200], addedinc[200], addedici[200];
	int minworkload[401];
	int wcet;
	k = 0;
	////////////////////////////////////////////////////////////////////////Ϊ���񸳳�ʼWCRT
	while (k<m) // WCRT = C for 'm' highest priority tasks
	{
		taskset[k][3] = taskset[k][0];
		k++;
	}
	while (k<n) // WCRT = D (worst-case) for all the rest of the tasks
	{
		taskset[k][3] = taskset[k][1];
		k++;
	}
	/////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////��ʼ������ɵ�����
	k = 0;
	// Initialize schedulability
	while (k<m) // UR_k_a = a, a = 1..C for 'm' highest priority tasks; set sched = 1
	{
		i = 1;
		while (i <= taskset[k][0])
		{
			wcrt[k][i] = i;
			i++;
		}
		schedulability[k] = 1;
		k++;
	}
	while (k<n) // for all other tasks, sched = 0
	{
		schedulability[k] = 0;
		k++;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////��ʼ��wcrt
	k = 0;
	// Initialize UR_k_0 = 0 for all n tasks
	while (k<n)
	{
		wcrt[k][0] = 0;
		k++;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////
	k = m;
	while (k<n)
	{
		//cout << "\tk = " << k << endl;
		//////////////////////////////////////////////////////////////////////////��ʼ��inc��oici Initialize inc and oici
		i = 0;
		while (i<n)
		{
			oinc[i] = 0;
			i++;
		}
		i = 0;
		while (i<n)
		{
			oici[i] = 0;
			i++;
		}
		i = 0;
		while (i < 401)
		{
			minworkload[i] = 0;
			i++;
		}
		//////////////////////////////////////////////////////////////////////////////////////////
		wcet = 1;
		//cout << "\t\twcrt:\t";
		//for(int k_index = 0; k_index < 8; ++k_index)
		//{
		//	cout << wcrt[k_index][0] << "\t";
		//}
		//cout << endl;
		/////////////////////////////////////////////////////////////////////����Ck=wcetʱk��WCRT
		while (wcet <= taskset[k][0]) // Theorem 5 - Calc WCRT when Ck = WCET
		{
			t = wcrt[k][wcet - 1] + 1;
			//if((curr_pr == 6) && (k == 6))
			//{
			//	cout << "\t\t\tk = " << k << "; wcet = " << wcet << "; t = " << t << endl;
			//}
			// CHANGE 1: no need to check if t <= D
			while (t >= 0)//<= taskset[k][1])
			{
				x = t - wcrt[k][wcet - 1];              //���䳤��; Expand the length = initial guess - UR_k_a
				i = 0;
				while (i<k)// Calc the interference of hp task i on task k
				{
					// Eq 2 - Within the interval length, cvalculate W^NC
					if (t%taskset[i][2]>taskset[i][0])	// t % Ti > Ci; so use Ci in Wnc formula
					{
						wnc[i] = int(t / taskset[i][2])*taskset[i][0] + taskset[i][0];
					}
					else 								// t % Ti < Ci; so use t % Ti in Wnc formular
					{
						wnc[i] = int(t / taskset[i][2])*taskset[i][0] + int(t%taskset[i][2]);
					}
					/////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////////////////////////////����t�����ڵ�W^{CI}
					if (t <= taskset[i][0])
					{
						wci[i] = t;
					}
					else
					{
						if ((t - taskset[i][0]) % taskset[i][2] > taskset[i][2] - taskset[i][3])
						{
							wci[i] = (int((t - taskset[i][0]) / taskset[i][2]) + 1)*taskset[i][0] + (t - taskset[i][0]) % taskset[i][2] - taskset[i][2] + taskset[i][3];
						}
						else
						{
							wci[i] = (int((t - taskset[i][0]) / taskset[i][2]) + 1)*taskset[i][0];
						}
					}
					/////////////////////////////////////////////////////////////////////////////////
					/////////////////////////////////////////////////////����WNC��WCI���������
					y = 0;
					maxaddedwnc = 0;
					while (y < taskset[i][2]&&y<=wcrt[k][wcet-1])
					{
						/////////////yʱ[0��wcrt[k][wcet-1])������i���nc����
						x1=(wcrt[k][wcet - 1] - y) / taskset[i][2] * taskset[i][0] + min(y, taskset[i][0]); ///////////////����i��[0��wcrt[k][wcet-1])�ϵ����nc���� 
						y1 = 0;///////////////////////////////////////////////����k��[wcrt[k][wcet-1]-y,wcrt[k][wcet-1])����Сִ��ʱ�䣻���y��Ӧ�Ĺ���������ֵ����,����i��[wcrt[k][wcet-1]-y,wcrt[k][wcet-1])��������k����С����
						while (wcrt[k][y1]+1 <= x + y)
						{
							y1++;
							if (y1 == wcet)
							{
								break;
							}
						}
						y1 = y1 - 1;
						if (y1 > y)
						{
							y++;
							continue;
						}
						if (y1 > taskset[i][0])
							y1 = taskset[i][0];
						if (y<= taskset[i][0])////////////////////////////////////////////////���x1������i��[0��wcrt[k][wcet-1])�϶�����k��������
						{
							x1 = x1 - y1;
						}
						else
						{
							y2 = 0;
							while (wcrt[k][y2] <= wcrt[k][wcet - 1] - y + taskset[i][0])
							{
								y2++;
							}
							y2 = y2 - 1;   
							x1 = x1-max(y1 + y2 - (wcet - 1), 0);
						}
						if (x1 < minworkload[i])
						{
							y++;
							continue;///////////////////////////////////wcrt[k][wcet-1]-y������Ϊ�ͷ�ʱ��
						}
						else
						{
							if (minworkload[i] <= (wcrt[k][wcet - 1] - y) / taskset[i][2] * taskset[i][0])/////////////y����Ӧ��Job����ִ��
							{
								if (taskset[i][3] <= y)
								{
									x2 = max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max((x + y - taskset[i][2]), 0) % taskset[i][2], taskset[i][0]);////x2ΪW�������
								}
								else
								{
									x2 = min(min(taskset[i][3]-y,taskset[i][0]-y1),x)+ max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max((x + y - taskset[i][2]), 0) % taskset[i][2], taskset[i][0]);
								}
							}
							else////////////////////////////////////////////////////[0,wcrt[k][wcet-1])��i����ִ��minworkload[i]-(wcrt[k][wcet - 1] - y) / taskset[i][2] * taskset[i][0]
							{
								y2 = minworkload[i] - (wcrt[k][wcet - 1] - y) / taskset[i][2] * taskset[i][0]+y1;////////////////////���y��Ӧ�Ĺ���������ֵ���أ���ô�˹�����[wcrt[k][wcet-1]-y,wcrt[k][wcet-1])�ϵ���С����
								if (taskset[i][3] <= y)
								{
									x2 = max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max((x + y - taskset[i][2]), 0) % taskset[i][2], taskset[i][0]);////x2ΪW�������
								}
								else
								{
									x2 = min(min(taskset[i][3] - y, max(taskset[i][0] - y2,0)), x) + max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max((x + y - taskset[i][2]), 0) % taskset[i][2], taskset[i][0]);
								}							
							}
							if (maxaddedwnc < x2)
							{
								maxaddedwnc = x2;
							}
						}
						y++;
					}
					y = 0;
					maxaddedwci = 0;
					while (y < taskset[i][2])
					{
						/////////////yʱ[0��wcrt[k][wcet-1])������i���nc����
						if (y <= wcrt[k][wcet - 1])
						{
							x1 = (wcrt[k][wcet - 1] - y) / taskset[i][2] * taskset[i][0] + min(y, taskset[i][0]) + max(min((wcrt[k][wcet - 1] - y) % taskset[i][2] - taskset[i][2] + taskset[i][3], taskset[i][0]), 0);
							y1 = 0;/////////////////////////    [ R_k^{wcet-1}-y, R_k^{wcet-1} )�ϵ���С����ʱ��  
							while (wcrt[k][y1]+1 <= x + y)
							{
								y1++;
								if (y1 == wcet)
								{
									break;
								}
							}
							y1 = y1 - 1;
							if (y1 > y)
							{
								y++;
								continue;
							}
							if (y1 > taskset[i][0])
								y1 = taskset[i][0];
							if (y <= taskset[i][0])
							{
								x1 = x1 - y1;
							}
							else
							{
								y2 = 0;
								while (wcrt[k][y2] <= wcrt[k][wcet - 1] - y + taskset[i][0])
								{
									y2++;
								}
								y2 = y2 - 1;
								x1 = x1 - max(y1 + y2 - (wcet - 1), 0);
							}
							if (x1 < minworkload[i])
							{
								y++;
								continue;///////////////////////////////////wcrt[k][wcet-1]-y������Ϊ�ͷ�ʱ��
							}
							else
							{
								if (minworkload[i] <= (wcrt[k][wcet - 1] - y) / taskset[i][2] * taskset[i][0] + max(min((wcrt[k][wcet - 1] - y) % taskset[i][2] - taskset[i][2] + taskset[i][3], taskset[i][0]), 0))/////////////y����Ӧ��Job����ִ��
								{
									if (taskset[i][3] <= y)
									{
										x2 = max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max(x + y - taskset[i][2], 0) % taskset[i][2], taskset[i][0]);
									}
									else
									{
										x2 = min(min(taskset[i][0]-y1, taskset[i][3]-y),x) + max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max(x + y - taskset[i][2], 0) % taskset[i][2], taskset[i][0]);
									}
								}
								else////////////////////////////////////////////////////[0,wcrt[k][wcet-1])��i����ִ��minworkload[i]-(wcrt[k][wcet - 1] - y) / taskset[i][2] * taskset[i][0]
								{
									y2 = minworkload[i] - (wcrt[k][wcet - 1] - y) / taskset[i][2] * taskset[i][0] - max(min((wcrt[k][wcet - 1] - y) % taskset[i][2] - taskset[i][2] + taskset[i][3], taskset[i][0]), 0)+y1;
									if (taskset[i][3] <= y)
									{
										x2 = max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max(x + y - taskset[i][2], 0) % taskset[i][2], taskset[i][0]);
									}
									else
									{
										x2 = min(min(max(taskset[i][0] - y2,0), taskset[i][3] - y),x) + max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max(x + y - taskset[i][2], 0) % taskset[i][2], taskset[i][0]);
									}
								}
								if (maxaddedwci < x2)
								{
									maxaddedwci = x2;
								}
							}
						}
						else///////////////////////y>wcrt[k][wcet-1]
						{
							y1 = wcet-1;///////////////////////////////////////////////////////////////�������i��wcrt[k][wcet-1]��ִ�У�����i��[0,wcrt[k][wcet-1])�ϵ���С����Ϊy1
							if (y1>taskset[i][0])
								y1 = taskset[i][0];
							if (taskset[i][3] >= y)
							{
								x1 = min(taskset[i][0], wcrt[k][wcet-1]);////////////////[0,wcrt[k][wcet-1])�ϵ������
								if (wcrt[k][wcet - 1] <= taskset[i][0])
									x1 = x1 - y1;////////////////////////////////////////[0,wcrt[k][wcet-1])�ϵ�������
								else
								{
									y2 = 0;
									while (wcrt[k][y2] <= taskset[i][0])
									{
										y2++;
										if (y2 == wcet)
										{
											break;
										}
									}
									y2 = y2 - 1;
									x1 = x1 - y2;
								}
							}
							else
							{
								x1 = max(min(taskset[i][0], taskset[i][3] + wcrt[k][wcet - 1] - y), 0);
								y2 = 0;
								while (wcrt[k][y2] <= x1)
								{
									y2++;
									if (y2 == wcet)
									{
										break;
									}
								}
								y2 = y2 - 1;
								x1 = x1 - y2;
							}
							if (x1 < minworkload[i])
							{
								y++;
								continue;///////////////////////////////////wcrt[k][wcet-1]-y������Ϊ�ͷ�ʱ��
							}
							else
							{
								y2 = minworkload[i] + y1;////////     �������i��wcrt[k][wcet-1]֮��ִ�У� ��Ҫ����+��С����
								if (taskset[i][3] <= y)
									x2 = max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max((x + y - taskset[i][2]), 0) % taskset[i][2], taskset[i][0]);        ////x2ΪW�������
								else
									x2 = min(min(max(min(taskset[i][0], taskset[i][3] - y + wcrt[k][wcet - 1]) - y2, 0),taskset[i][3]-y), x)+max((x + y - taskset[i][2]), 0) / taskset[i][2] * taskset[i][0] + min(max((x + y - taskset[i][2]), 0) % taskset[i][2], taskset[i][0]);
								if (maxaddedwci < x2)
								{
									maxaddedwci = x2;
								}
							}
						}
						y++;
					}
					/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////�������W=I1+W2
					if (wnc[i]>oinc[i] + maxaddedwnc)
						wnc[i] = oinc[i] + maxaddedwnc;
					if (wci[i] > oici[i] + maxaddedwci)
						wci[i] = oici[i] + maxaddedwci;
					if (wnc[i] < oinc[i] + maxaddedwnc)
						wnc[i] = wnc[i] ;
					if (wci[i] < oici[i] + maxaddedwci)
						wnc[i] = wnc[i] ;
					/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					////////////////////////////////////////////////////////////////����t�����ڵ�I
					// if((curr_pr == 6) && (k == 6))
					// {
					// 	cout << "\t\tLINE 2682 t = " << t << endl;
					// }
					ce = t - wcet + 1;
					if (wci[i] > t - wcet + 1)
					{
						ici[i] = t - wcet + 1;
				//		limitednumber = limitednumber + 1;
					}
					else
					{
						ici[i] = wci[i];
					}
					if (wnc[i] > t - wcet + 1)
					{
						inc[i] = t - wcet + 1;
				//		limitednumber = limitednumber + 1;
					}
					else
					{
						inc[i] = wnc[i];
					}
					//////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////////////////////////////��¼i��t�ϵ�������
					//////////////////////////////////////////////////////////////////////////////
					diffw[i] = ici[i] - inc[i];///////////////////////////////////////����I^{diff}
					i++;
				}////////////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////�����ܸ���
				i = 0;
				interference = 0;
				while (i<k)
				{
					interference = interference + inc[i];
					i++;
				}
				j = 0;
				while (j<m - 1)
				{
					l = 0;
					o = 0;
					while (l<k)
					{
						if (diffw[o]<diffw[l])
						{
							o = l;
						}
						l++;
					}
					interference = interference + diffw[o];
					diffw[o] = 0;
					j++;
				}
				////////////////////////////////////////////////////////////////////////////////
				interference = int(interference / m);
				//cout << "\t\tinterference = " << interference << "; wcet = " << wcet << "; interference + wcet = " << interference + wcet << "; t = " << t << endl;
				if (interference + wcet <= t)
				{
					wcrt[k][wcet] = t;
					/////////////////////////////////////////////////////����ÿ��������t�ϵ���С����
					j = 0;
					while (j < k)
					{
						if (ici[j]>interference)
						{
							ici[j] = interference;
						}
						if (inc[j]>interference)
							inc[j] = interference;
						j++;
					}
					i = 0;
					while (i < k)
					{
						j = 0;
						sumi= 0;
						while (j<k)
						{
							if (j != i)
							{
								sumi = sumi+ inc[i];
							}							
							j++;
						}
						j = 0;
						while (j < k)
						{
							diffi[j] = max(ici[j] - inc[j],0);
							j++;
						}
						j = 0;
						while (j<m - 1)
						{
							l = 0;
							o = 0;
							while (l<k)
							{
								if (l != i)
								{
									if (diffw[o] < diffw[l])
									{
										o = l;
									}
								}
								l++;
							}
							sumi = sumi + diffi[o];
							diffi[o] = 0;
							j++;
						}
						minworkload[i] = MY_MAX(interference*m - sumi,0);
						i++;
					}
					/////////////////////////////////////////////////////////////////////////////////
					// a = C; done computing UR_k, RT = t
					if (wcet == taskset[k][0])
					{
						schedulability[k] = 1;
						taskset[k][3] = t;
					}
					// cout << "\t\twcrt after a = Ck; k = " << k << "; a = " << wcet << "; t = " << t << ":\t";
					// for(int k_index = 0; k_index < 8; ++k_index)
					// {
					// 	cout << wcrt[k_index][0] << "\t";
					// }
					// cout << endl;

					break;
				}
				else
				{
					t = interference + wcet;
				}
				// cout << "\t\twcrt after inner loop:\t";
				// for(int k_index = 0; k_index < 8; ++k_index)
				// {
				// 	cout << wcrt[k_index][0] << "\t";
				// }
				// cout << endl;
			}
			// CHANGE 2: continue computing R even if R > D
			// if (t>taskset[k][1])  ////////////////////Rk>Dk
			// {
			// 	wcrt[k][wcet] = taskset[k][1];
			// 	taskset[k][3] = taskset[k][1];
			// 	schedulability[k] = 0;
			// 	break;
			// }
			i = 0;
			while (i < k)
			{
				if (inc[i] <= t - wcet)
					oinc[i] = inc[i];
				else
					oinc[i] = t - wcet;
				if (oici[i] <= t - wcet)
					oici[i] = ici[i];
				else
					oici[i] = t-wcet;
				i++;
			}
			//cout << "k = " << k << " t = " << t << " wcet = " << wcet << "\n";
			wcet++;
		}
		////////////////////////////////////////////////////////////////////////////////////////CHANGE 3
		//cout << "Outside: k = " << k << " t = " << t << " wcet = " << wcet << "\n";
		//CHANGE 3: stop and return response time when priority reaches curr_pr
		// if (schedulability[k] == 0)
		// {
		// 	break;
		// }
		if(k == curr_pr)
		{
			//cout << "a = " << wcet << "; response time: " << t << endl;
			//cout << "This is executed; k = " << k << "; curr_pr = " << curr_pr << "; t = " << t << endl; 
			return t; 
		}
		////////////////////////////////////////////////////////////////////////////////////////
		k++;
	}
	//////////////////////////////////////////////////////////////////////////////////�жϿɵ�����
	cout << "Control never reaches here\n";
	i = 0;
	result = 1;
	while (i<n)
	{
		if (schedulability[i] == 0)
		{
			result = 0;
			break;
		}
		i++;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////
	return result;
}

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