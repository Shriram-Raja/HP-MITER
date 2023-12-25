#pragma once
#include "MultiProcessorSystem.h"
#include "SchedulabilityAnalysis.h"
#include "GeneralUnschedCoreComputer.hpp"
#ifdef MILP
#include "GeneralFocusRefinement.h"
#endif
#include <thread>
#include <mutex>
#include "BitMask.h"
namespace MPScheduling {	
	class MPSchedulingUnschedCoreComputer : public GeneralUnschedCoreComputer_MonoInt<MultiProcessorSystem, pair<int, int>, ValueType> {
	public:
		struct RelaxedSchedulablePAList {
			mutex m_cAccessLock;
			deque<PriorityAssignment> m_dequePAs;
			void push(PriorityAssignment & rcPA);
			vector<PriorityAssignment> pop();
		};
	protected:
		RelaxedSchedulablePAList m_cRelaxedSchedulablePAList;
		typedef unordered_map<int, BitMaskBasedSubSupsetQuery> SchedCache;
		unordered_map<int, BitMaskBasedSubSupsetQuery> m_cSchedCachedForMUSConversion;
		unordered_map<int, BitMaskBasedSubSupsetQuery> m_cUnschedCachedForMUSConversion;
		vector<int> m_vectorOrder;
	public:
		MPSchedulingUnschedCoreComputer();
		MPSchedulingUnschedCoreComputer(SystemType & rcSystem);
		~MPSchedulingUnschedCoreComputer() {}
		virtual SchedCondContent SchedCondUpperBound(const SchedCondType & rcSchedCondType);
		virtual SchedCondContent SchedCondLowerBound(const SchedCondType & rcSchedCondType);
		virtual Monotonicity SchedCondMonotonicity(const SchedCondType & rcSchedCondType);
		virtual bool IsSchedulable(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed);
		virtual bool IsSchedulable(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed, PriorityAssignment & rcPriorityAssignment);
		RelaxedSchedulablePAList & getRelaxedSchedulablePAList() { return m_cRelaxedSchedulablePAList; }
	protected:		
		void ExtractRTEAssginment(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed, vector<vector<ValueType>> & rvectorRTERange);		
		int ConvertToMUS(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed, SchedConds & rcCore);		
		int isSchedulable(vector<vector<ValueType>> & rvectorRTERange, PriorityAssignment & rcPA, int iStartLevel);
		int isSchedulableWithSchedCache(vector<vector<ValueType>> & rvectorRTERange, PriorityAssignment & rcPA, int iStartLevel, SchedCache & rcSchedCache);
		int ConvertToMUS(vector<vector<ValueType>> & rvectorRTERange);		
		void ComputeMaxUnschedRTEUB(int iTaskIndex, vector<vector<ValueType>> & rvectorRTERange, pair<PriorityAssignment, int> & rcPAInfo);
		void ComputeMaxUnschedRTELB(int iTaskIndex, vector<vector<ValueType>> & rvectorRTERange, pair<PriorityAssignment, int> & rcPAInfo);
		void ComputeMaxUnschedRTELB_SchedCache(int iTaskIndex, vector<vector<ValueType>> & rvectorRTERange, pair<PriorityAssignment, int> & rcPAInfo);
		int ComputeMaxUnschedRTELB_isSchedulable(int iTaskIndex, vector<vector<ValueType>> & rvectorRTERange, PriorityAssignment & rcPA, int iStartLevel, vector<int> & rvectorOrder, int iWorkType);
		int ComputeMaxUnschedRTELB_isSchedulable_SchedCache(int iTaskIndex, vector<vector<ValueType>> & rvectorRTERange, PriorityAssignment & rcPA, int iStartLevel, vector<int> & rvectorOrder, int iWorkType, SchedCache & rcSchedCache);
		int ComputeMaxUnschedRTELB_isSchedulable_SchedCache(int iTaskIndex, vector<vector<ValueType>> & rvectorRTERange, PriorityAssignment & rcPA, int iStartLevel, vector<int> & rvectorOrder, int iWorkType, SchedCache & rcSchedCache, SchedCache & rcUnschedCache, SchedCache & rcNewUnschedCache);
		void AddSchedCached(PriorityAssignment & rcPA, int iLB, SchedCache & rcSchedCachedForMUSConversion);		
	};
	
#ifdef MILP	
	class MPSchedulingFR : public GeneralFocusRefinement_MonoInt<MPSchedulingUnschedCoreComputer> {
	public:
		struct Configuration {
			int m_iMaximizeMinimumSlackObjective;
			int m_iUseAsyncTest;
			Configuration();
		};
	protected:
		vector<thread> m_vectorAsynchronousTestThread;
		struct AsynchronousTestControl {
			deque<PriorityAssignment> m_dequeFound;
			mutex m_cLock;
			bool m_bTerminate;
			bool empty();
			AsynchronousTestControl() :m_bTerminate(false) {}
			void push(PriorityAssignment & rcPA);
			void terminate();
			bool isTerminate();
			PriorityAssignment top();
		} m_cAsyncTestControl;		
		Configuration m_cConfiguration;
	public:
		MPSchedulingFR() {}
		MPSchedulingFR(SystemType & rcSystem);
		int FocusRefinement(int iDisplay, double dTimeout = 1e74);
		Configuration & getConfiguration() { return m_cConfiguration; }
	protected:
		void CreateSchedCondTypeVars(const IloEnv & rcEnv, SchedCondTypeVars & rcVars);
		void GenerateLBUBEqualConst(const IloEnv & rcEnv, SchedCondTypeVars & rcVars, IloRangeArray & rcConst);
		IloObjective getMaximizeMinimumSlackObjective(const IloEnv & rcEnv, SchedCondTypeVars & rcVars, IloRangeArray & rcConst);
		int SolveFR(int iDisplay, SchedConds & rcSchedCondFixed, UnschedCores & rcGenUC, double dObjLB, double dObjUB, double dTimeout = 1e74);		
		void StartAsynchronousTest(int iThreadNum);		
		void CloseAsynchronousTest();
	};
	void DebugFR();
#endif

	void TestMPSchedulingFR();
	void TestOptimizedUnschedCoreComputer();

	class ExhaustiveSearch_GuanRTA {
	public:
		unsigned long long count;
		ExhaustiveSearch_GuanRTA() {}
		pair<bool, PriorityAssignment> operator () (MultiProcessorSystem & rcSystem);
	protected:
		bool recur(int iLevel, MultiProcessorSystem & rcSystem, PriorityAssignment & rcPA);		
	};	
}