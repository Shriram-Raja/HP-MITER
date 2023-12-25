#pragma once
#include "MultiProcessorSystem.h"

// OPA-Compatible "tests"
#define DA		0
#define DA_LC	1
#define C_RTA	2

namespace MPScheduling {
	class ResponseTimeAnalysis_Guan {
	public:
	protected:
	public:
		ResponseTimeAnalysis_Guan();
		~ResponseTimeAnalysis_Guan() {}
		bool ComputeAllRTGivenPA(MultiProcessorSystem & rcSystem, const PriorityAssignment & rcPA, vector<ValueType> & rvectorResponseTime);
		pair<bool, vector<ValueType>> ComputeAllRTGivenPA(MultiProcessorSystem & rcSystem, const PriorityAssignment & rcPA);
		ValueType RT(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA);
		ValueType RTTerminateOnBound(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA, ValueType cBound);
		ValueType RT(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA, ValueType cLB, ValueType cUB);
		ValueType RTNoCI(MultiProcessorSystem & rcSystem, int iTaskIndex, const PriorityAssignment & rcPA);	
		bool CRTA(MultiProcessorSystem & rcSystem, int kTaskIndex, PriorityAssignment & rcPA);
	protected:
		ValueType CarryinInterference(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, ValueType dResponseTime);
		ValueType WCI(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, ValueType dResponseTime);
		ValueType WNC(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t);		
		ValueType RHS(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA);						
	};

	// *********************** Aug.5 Patch one: Add ZLL class for experiments ******************************
	class ResponseTimeAnalysis_ZLL {
	public:
			ResponseTimeAnalysis_ZLL();
			~ResponseTimeAnalysis_ZLL() {};
			// ----------------------------------- Part 1: GSYY related computation -----------------------------------------
			bool ComputeAllRTGivenPA(MultiProcessorSystem & rcSystem, const PriorityAssignment & rcPA, vector<ValueType> & rvectorResponseTime);
			pair<bool, vector<ValueType>> ComputeAllRTGivenPA(MultiProcessorSystem & rcSystem, const PriorityAssignment & rcPA);
			ValueType RT(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA);
			
			ValueType RT(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA, ValueType cLB, ValueType cUB);
			ValueType RTNoCI(MultiProcessorSystem & rcSystem, int iTaskIndex, const PriorityAssignment & rcPA);		
	protected:
			ValueType CarryinInterference(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, ValueType dResponseTime);
			ValueType WCI(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, ValueType dResponseTime);
			ValueType WNC(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t);		
			ValueType RHS(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA);
			ValueType ComputeOmega(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA);
			// ----------------------------------- End Part 1: GSYY related computation -----------------------------------------

	public:
		//  Part 2: ZLL related computation, Compute carry-in workload beta based on Lemma 1, Lemma 2 and Eqn.(20)
			ValueType ComputeMIN(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA); // Lemma 1
			ValueType ComputeMAX(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA); // Lemma 2
			ValueType ComputeBeta(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA); // Eqn.(20)
			
			ValueType WCI_ZLL(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA); // Eqn.(20)
			ValueType RHS_ZLL(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA);
			ValueType RTTerminateOnBound(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<ValueType> & rvectorResponseTime, const PriorityAssignment & rcPA, ValueType cBound);
	};
	// *********************** Aug.5 Patch one: Add ZLL class for experiments ******************************

	/**
	 * @brief Implementation of EPE (EMSOFT 2021) schedulability analysis
	 * @date  Mar. 30, 2022
	 * @version v1 for initial test 
	 */

	class ResponseTimeAnalysis_EPE{
	public: 
			ResponseTimeAnalysis_EPE();
			~ResponseTimeAnalysis_EPE() {};

	public:
		bool GSYY_ComputeAllRTGivenPA(MultiProcessorSystem & rcSystem, const PriorityAssignment & rcPA, vector<vector<ValueType>> & rvectorResponseTime);
		pair<bool, vector<vector<ValueType>>> GSYY_ComputeAllRTGivenPA(MultiProcessorSystem & rcSystem, const PriorityAssignment & rcPA);
		ValueType GSYY_RT(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<vector<ValueType>> & rvectorResponseTime, const PriorityAssignment & rcPA);
		ValueType GSYY_RTTerminateOnBound(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<vector<ValueType>> & rvectorResponseTime, const PriorityAssignment & rcPA, ValueType cBound);
		ValueType GSYY_RT(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<vector<ValueType>> & rvectorResponseTime, const PriorityAssignment & rcPA, ValueType cLB, ValueType cUB);
		ValueType GSYY_RTNoCI(MultiProcessorSystem & rcSystem, int iTaskIndex, const PriorityAssignment & rcPA);		

		// GSYY Overloaded
		ValueType GSYY_RTTerminateOnBound(MultiProcessorSystem & rcSystem, int iTaskIndex, vector<vector<ValueType>> & rvectorResponseTime, const PriorityAssignment & rcPA, ValueType cBound, int a);

		// GSYY2 Interfaces
		pair<bool, vector<ValueType>> ResponseTimeAnalysis_EPE::ComputeRTAGSYY2(MultiProcessorSystem & rcSystem, const PriorityAssignment & rcPA);

		// GSYY2 (EPE) interface to compute response time of task at priority 'curr_pr'
		ValueType ResponseTimeAnalysis_EPE::ComputeRTAGSYY2(MultiProcessorSystem & rcSystem, const PriorityAssignment & rcPA, int curr_pr);

	protected:
		// ----- Part 1: GSYY base -----
		ValueType GSYY_CarryinInterference(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, ValueType dResponseTime);
		ValueType GSYY_WCI(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, ValueType dResponseTime);
		ValueType GSYY_WNC(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t);		
		ValueType GSYY_RHS(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, vector<vector<ValueType>> & rvectorResponseTime, const PriorityAssignment & rcPA);		

		// ----- Part 2a: GSYY Overloaded -----
		ValueType GSYY_RHS(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t, vector<vector<ValueType>> & rvectorResponseTime, const PriorityAssignment & rcPA, int a);
	};



	class DeadlineAnalysis {
	public:

	protected:

	public:
		DeadlineAnalysis() {}
		~DeadlineAnalysis() {}
		bool isSchedulable(MultiProcessorSystem & rcSystem, int iTaskIndex, PriorityAssignment & rcPA);
		bool isSchedulable_LC(MultiProcessorSystem & rcSystem, int iTaskIndex, PriorityAssignment & rcPA); // DA-LC
		bool Audsley(MultiProcessorSystem & rcSystem, PriorityAssignment & rcPA, int sched_test = 0);
		ValueType Interference(MultiProcessorSystem & rcSystem, int iTaskIndex, PriorityAssignment & rcPA);
	protected:		
		ValueType NDi(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t);
		ValueType WDi(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t);
		ValueType IDi(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType Dk, ValueType Ck);

		ValueType WNCi(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType t);
		ValueType INCi(MultiProcessorSystem & rcSystem, int iTaskIndex, ValueType Dk, ValueType Ck);		
	};
}
