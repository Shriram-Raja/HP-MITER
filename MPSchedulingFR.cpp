#include "stdafx.h"
#include "MPSchedulingFR.h"
#include <algorithm>
#include "SchedulabilityAnalysis.h"
#include "MultiProcessorSystem.h"
#include<list>


namespace MPScheduling {
	ostream & operator << (ostream & os, const MPSchedulingUnschedCoreComputer::SchedCondType & rcType) {
		os << '{' << rcType.first << ", " << rcType.second << '}';
		return os;
	}

	istream & operator >> (istream & is, MPSchedulingUnschedCoreComputer::SchedCondType & rcType) {
		while (is.peek() == ' ') is.get();
		char ch = is.get();
		assert(ch == '{');
		is >> rcType.first;
		assert(is.get() == ',');
		assert(is.get() == ' ');
		is >> rcType.second;
		assert(is.get() == '}');
		return is;
	}

	//RelaxedSchedulablePAList
	void MPSchedulingUnschedCoreComputer::RelaxedSchedulablePAList::push(PriorityAssignment & rcPA) {
		m_cAccessLock.lock();
		m_dequePAs.push_back(rcPA);
		m_cAccessLock.unlock();
	}

	vector<PriorityAssignment> MPSchedulingUnschedCoreComputer::RelaxedSchedulablePAList::pop() {
		m_cAccessLock.lock();
		vector<PriorityAssignment> ret;
		if (m_dequePAs.empty() == false) {
			ret.push_back(m_dequePAs.back());
			m_dequePAs.pop_back();
		}
		m_cAccessLock.unlock();
		return ret;
	}

	MPSchedulingUnschedCoreComputer::MPSchedulingUnschedCoreComputer() {		
	}

	MPSchedulingUnschedCoreComputer::MPSchedulingUnschedCoreComputer(SystemType & rcSystem)
		:GeneralUnschedCoreComputer_MonoInt(rcSystem) {
		m_vectorOrder = vector<int>(m_pcTaskSet->getNumTasks());
		for (int i = 0; i < m_pcTaskSet->getNumTasks(); i++) m_vectorOrder[i] = i;
		sort(m_vectorOrder.begin(), m_vectorOrder.end(), [&](int lhs, int rhs)->bool {return rcSystem.getTask(lhs).D() > rcSystem.getTask(rhs).D();});
	}

	MPSchedulingUnschedCoreComputer::SchedCondContent MPSchedulingUnschedCoreComputer::SchedCondUpperBound(const SchedCondType & rcSchedCondType) {
		return min(m_pcTaskSet->getTask(rcSchedCondType.first).D(), m_pcTaskSet->getTask(rcSchedCondType.first).T());
	}

	MPSchedulingUnschedCoreComputer::SchedCondContent MPSchedulingUnschedCoreComputer::SchedCondLowerBound(const SchedCondType & rcSchedCondType) {
		return m_pcTaskSet->getTask(rcSchedCondType.first).C();
	}

	MPSchedulingUnschedCoreComputer::Monotonicity MPSchedulingUnschedCoreComputer::SchedCondMonotonicity(const SchedCondType & rcSchedCondType) {
		return rcSchedCondType.second == 1 ? Monotonicity::Up : Monotonicity::Down;
	}

	void MPSchedulingUnschedCoreComputer::ExtractRTEAssginment(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed, vector<vector<ValueType>> & rvectorRTERange) {
		if (rvectorRTERange.size() != 2) rvectorRTERange = vector<vector<ValueType>>(2, vector<ValueType>(m_pcTaskSet->getNumTasks()));
		for (int i = 0; i < m_pcTaskSet->getNumTasks(); i++) {
			rvectorRTERange[0][i] = m_pcTaskSet->getTask(i).C();
			rvectorRTERange[1][i] = m_pcTaskSet->getTask(i).D();
		}
		auto update = [&](SchedCond ele) {
			rvectorRTERange[ele.first.second][ele.first.first] = (ele.first.second == 0) ?
				max(rvectorRTERange[ele.first.second][ele.first.first], ele.second) : min(rvectorRTERange[ele.first.second][ele.first.first], ele.second);
		};
		for (auto ele : rcSchedCondsFlex) update(ele);
		for (auto ele : rcSchedCondsFixed) update(ele);
	}

	bool MPSchedulingUnschedCoreComputer::IsSchedulable(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed) {
		PriorityAssignment cDummyPA(*m_pcTaskSet);
		return IsSchedulable(rcSchedCondsFlex, rcSchedCondsFixed, cDummyPA);
	}
	double isSchedulableTime = 0;
	bool MPSchedulingUnschedCoreComputer::IsSchedulable(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed, PriorityAssignment & rcPriorityAssignment) {
		vector<vector<ValueType>> vectorRTERange;
		int iTaskNum = m_pcTaskSet->getNumTasks();
		ExtractRTEAssginment(rcSchedCondsFlex, rcSchedCondsFixed, vectorRTERange);
		int iPriority = iTaskNum - 1;
		for (; iPriority >= 0; iPriority--)
		{
			bool bAssigned = false;
			if (rcPriorityAssignment.getTaskByPriority(iPriority) != -1) continue;
			for (int i = 0; i < iTaskNum && bAssigned == false; i++)
			{
				if ((rcPriorityAssignment.getPriorityByTask(i) != -1)) continue;
				rcPriorityAssignment.setPriority(i, iPriority);
				auto rt = ResponseTimeAnalysis_Guan().RTTerminateOnBound(*m_pcTaskSet, i, vectorRTERange[0], rcPriorityAssignment, vectorRTERange[1][i]);
				bAssigned = rt <= vectorRTERange[1][i];
				if (!bAssigned) rcPriorityAssignment.unset(i);
			}
			if (!bAssigned) return false;
		}

		return true;
	}

	int MPSchedulingUnschedCoreComputer::ConvertToMUS(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed, SchedConds & rcCore) {
		//return GeneralUnschedCoreComputer_MonoInt::ConvertToMUS(rcSchedCondsFlex, rcSchedCondsFixed, rcCore);
		vector<vector<ValueType>> vectorRTERange;
		ExtractRTEAssginment(rcSchedCondsFlex, rcSchedCondsFixed, vectorRTERange);		
		if (ConvertToMUS(vectorRTERange) == 0) return 0;
		for (int i = 0; i < m_pcTaskSet->getNumTasks(); i++) {
			if (vectorRTERange[0][i] != m_pcTaskSet->getTask(i).C()) rcCore[SchedCondType(i, 0)] = vectorRTERange[0][i];
			if (vectorRTERange[1][i] != m_pcTaskSet->getTask(i).D()) rcCore[SchedCondType(i, 1)] = vectorRTERange[1][i];
		}
		return 1;
	}
	
	int MPSchedulingUnschedCoreComputer::isSchedulable(vector<vector<ValueType>> & rvectorRTERange, PriorityAssignment & rcPA, int iStartLevel) {
		return isSchedulableWithSchedCache(rvectorRTERange, rcPA, iStartLevel, m_cSchedCachedForMUSConversion);
		int iTaskNum = m_pcTaskSet->getNumTasks();
		int iPriority = iStartLevel;
		for (; iPriority >= 0; iPriority--)
		{
			bool bAssigned = false;
			for (int i = 0; i < iTaskNum && bAssigned == false; i++)
			{
				if ((rcPA.getPriorityByTask(i) > iPriority)) continue;
				rcPA.unset(i);
				rcPA.setPriority(i, iPriority);
				auto rt = ResponseTimeAnalysis_Guan().RTTerminateOnBound(*m_pcTaskSet, i, rvectorRTERange[0], rcPA, rvectorRTERange[1][i]);
				bAssigned = rt <= rvectorRTERange[1][i];
				if (!bAssigned) rcPA.unset(i);
			}
			if (!bAssigned) break;
		}
		if (iPriority == -1) m_cRelaxedSchedulablePAList.push(rcPA);
		return iPriority;
	}

	int MPSchedulingUnschedCoreComputer::isSchedulableWithSchedCache(vector<vector<ValueType>> & rvectorRTERange, PriorityAssignment & rcPA, int iStartLevel, SchedCache & rcSchedCache) {		
		int iTaskNum = m_pcTaskSet->getNumTasks();
		int iPriority = iStartLevel;
		BitMask cLPSet(iTaskNum);
		for (int i = iTaskNum - 1; i > iPriority; i--) cLPSet.set(rcPA.getTask(i));		
		for (; iPriority >= 0; iPriority--)
		{
			bool bAssigned = false;
			for (int i = 0; i < iTaskNum && bAssigned == false; i++) {
				if ((rcPA.getPriorityByTask(i) > iPriority)) continue;
				if (bAssigned = (rcSchedCache[i].existSubset(cLPSet).first)) {
					rcPA.setPriority(i, iPriority); 
					cLPSet.set(i);
					break;
				}
			}
			for (int i = 0; i < iTaskNum && bAssigned == false; i++)
			{
				if ((rcPA.getPriorityByTask(i) > iPriority)) continue;
				rcPA.unset(i);
				rcPA.setPriority(i, iPriority);
				auto rt = ResponseTimeAnalysis_Guan().RTTerminateOnBound(*m_pcTaskSet, i, rvectorRTERange[0], rcPA, rvectorRTERange[1][i]);				
				bAssigned = rt <= rvectorRTERange[1][i];
				if (!bAssigned) rcPA.unset(i);
				else cLPSet.set(i);
			}
			if (!bAssigned) break;			
		}
		if (iPriority == -1) m_cRelaxedSchedulablePAList.push(rcPA);
		return iPriority;
	}

	int MPSchedulingUnschedCoreComputer::ConvertToMUS(vector<vector<ValueType>> & rvectorRTERange) {
		Timer cTimer; cTimer.Start();
		int iTaskNum = m_pcTaskSet->getNumTasks();
		m_cSchedCachedForMUSConversion.clear();
		pair<PriorityAssignment, int> cPAInfo(PriorityAssignment(*m_pcTaskSet), 0);
		cPAInfo.second = isSchedulable(rvectorRTERange, cPAInfo.first, iTaskNum - 1);
		if (cPAInfo.second == -1) return 0;
		for (int i = 0; i < iTaskNum; i++) {
			//ComputeMaxUnschedRTELB(i, rvectorRTERange, cPAInfo);
			ComputeMaxUnschedRTELB_SchedCache(i, rvectorRTERange, cPAInfo);
			ComputeMaxUnschedRTEUB(i, rvectorRTERange, cPAInfo);
		}
		cTimer.Stop();
		isSchedulableTime += cTimer.getElapsedTime_ms();
		return 1;
	}

	void MPSchedulingUnschedCoreComputer::ComputeMaxUnschedRTEUB(int iTaskIndex, vector<vector<ValueType>> & rvectorRTERange, pair<PriorityAssignment, int> & rcPAInfo) {
		PriorityAssignment & rcPA = rcPAInfo.first;
		if (rcPA.getPriority(iTaskIndex) != -1) {
			rvectorRTERange[1][iTaskIndex] = m_pcTaskSet->getTask(iTaskIndex).D();
			return;
		}

		rcPA.setPriority(iTaskIndex, rcPAInfo.second);
		ValueType cRT = ResponseTimeAnalysis_Guan().RTTerminateOnBound(*m_pcTaskSet, iTaskIndex, rvectorRTERange[0], rcPA, m_pcTaskSet->getTask(iTaskIndex).D());
		if (cRT > m_pcTaskSet->getTask(iTaskIndex).D()) {
			rcPA.unset(iTaskIndex);
			rvectorRTERange[1][iTaskIndex] = m_pcTaskSet->getTask(iTaskIndex).D();
			return;
		}

		int iStatus = isSchedulable(rvectorRTERange, rcPA, rcPAInfo.second - 1);
		if (iStatus == -1) {
			rvectorRTERange[1][iTaskIndex] = cRT - 1;
			AddSchedCached(rcPA, rcPAInfo.second - 1, m_cSchedCachedForMUSConversion);
			for (int i = rcPAInfo.second; i >= 0; i--) rcPA.unset(rcPA.getTask(i));
		}
		else {
			rvectorRTERange[1][iTaskIndex] = m_pcTaskSet->getTask(iTaskIndex).D();
			rcPAInfo.second = iStatus;
		}
		return;
	}

	double timeComputeMaxUnschedRTELB = 0;
	int numComputeMaxUnschedRTELB = 0;
	void MPSchedulingUnschedCoreComputer::ComputeMaxUnschedRTELB(int iTaskIndex, vector<vector<ValueType>> & rvectorRTERange, pair<PriorityAssignment, int> & rcPAInfo) {
		Timer cTimer; cTimer.Start();
		int iTaskNum = m_pcTaskSet->getNumTasks();
		PriorityAssignment & rcPA = rcPAInfo.first;
		ValueType cUnschedBound = rvectorRTERange[0][iTaskIndex], cSchedBound = m_pcTaskSet->getTask(iTaskIndex).C();
		rvectorRTERange[0][iTaskIndex] = m_pcTaskSet->getTask(iTaskIndex).C();
		if (rcPA.getPriority(iTaskIndex) != -1) {
			cTimer.Stop(); timeComputeMaxUnschedRTELB += cTimer.getElapsedTime_ms();
			return;
		}

		vector<int> vectorOrder = { iTaskIndex };
		for (int i = 0; i < iTaskNum; i++) if (rcPA.getPriority(i) <= rcPAInfo.second && i != iTaskIndex) vectorOrder.push_back(i);
		sort(vectorOrder.begin() + 1, vectorOrder.end(), [&](int lhs, int rhs)->bool {return m_pcTaskSet->getTask(lhs).D() > m_pcTaskSet->getTask(rhs).D();});
		assert(vectorOrder.size() == rcPAInfo.second + 1);

		int iStatus = ComputeMaxUnschedRTELB_isSchedulable(iTaskIndex, rvectorRTERange, rcPA, rcPAInfo.second, vectorOrder, 0);
		if (iStatus != -1) {
			cTimer.Stop(); timeComputeMaxUnschedRTELB += cTimer.getElapsedTime_ms();
			rcPAInfo.second = iStatus;
			return;
		}

		auto cExtractOrder = [&]()->void {
			for (int i = rcPAInfo.second; i >= 0; i--) {
				assert(rcPA.getTask(i) != -1);
				vectorOrder.at(rcPAInfo.second - i) = rcPA.getTask(i);
			}
		};
		cExtractOrder();
		AddSchedCached(rcPA, rcPA.getPriority(iTaskIndex), m_cSchedCachedForMUSConversion);

		PriorityAssignment cUnschedPA = rcPAInfo.first;
		int iUnschedFrontier = rcPAInfo.second;
		while (cUnschedBound - cSchedBound > 1) {
			rvectorRTERange[0][iTaskIndex] = (cSchedBound + cUnschedBound) >> 1;
			int iStatus = ComputeMaxUnschedRTELB_isSchedulable(iTaskIndex, rvectorRTERange, rcPA, rcPAInfo.second, vectorOrder, 1);
			if (iStatus == -1) {
				cSchedBound = rvectorRTERange[0][iTaskIndex];
				cExtractOrder();
				AddSchedCached(rcPA, rcPA.getPriority(iTaskIndex), m_cSchedCachedForMUSConversion);
			}
			else {
				cUnschedBound = rvectorRTERange[0][iTaskIndex];
				cUnschedPA = rcPA;
				iUnschedFrontier = iStatus;
			}
		}
		assert(iUnschedFrontier != -1);
		rvectorRTERange[0][iTaskIndex] = cUnschedBound;
		for (int i = iUnschedFrontier; i >= 0; i--) cUnschedPA.unset(cUnschedPA.getTask(i));
		rcPAInfo = { cUnschedPA, iUnschedFrontier };
		cTimer.Stop();
		timeComputeMaxUnschedRTELB += cTimer.getElapsedTime_ms();
	}

	int MPSchedulingUnschedCoreComputer::ComputeMaxUnschedRTELB_isSchedulable(int iTaskIndex, vector<vector<ValueType>> & rvectorRTERange, PriorityAssignment & rcPA, int iStartLevel, vector<int> & rvectorOrder, int iWorkType) {
		return ComputeMaxUnschedRTELB_isSchedulable_SchedCache(iTaskIndex, rvectorRTERange, rcPA, iStartLevel, rvectorOrder, iWorkType, m_cSchedCachedForMUSConversion);
		int iTaskNum = m_pcTaskSet->getNumTasks();
		list<int> listCandidates(rvectorOrder.begin(), rvectorOrder.end());
		for (auto cand : rvectorOrder) rcPA.unset(cand);
		int iPriority = iStartLevel;
		for (; iPriority >= 0; iPriority--)
		{
			bool bAssigned = false;
			auto iter = listCandidates.begin();
			for (; iter != listCandidates.end(); iter++) {
				int iCand = *iter;
				rcPA.setPriority(iCand, iPriority);
				auto rt = ResponseTimeAnalysis_Guan().RTTerminateOnBound(*m_pcTaskSet, iCand, rvectorRTERange[0], rcPA, rvectorRTERange[1][iCand]);
				bAssigned = rt <= rvectorRTERange[1][iCand];
				if (!bAssigned) rcPA.unset(iCand);
				else break;
			}
			if (!bAssigned) break;
			listCandidates.erase(iter);

			if (iWorkType == 1 && rcPA.getTask(iPriority) == iTaskIndex) {
				for (auto remaining : listCandidates) rcPA.setPriority(remaining, --iPriority);
				assert(rcPA.getSize() == iTaskNum);
				return -1;
			}
		}
		if (iPriority == -1) m_cRelaxedSchedulablePAList.push(rcPA);
		return iPriority;
	}

	int MPSchedulingUnschedCoreComputer::ComputeMaxUnschedRTELB_isSchedulable_SchedCache(int iTaskIndex, vector<vector<ValueType>> & rvectorRTERange, PriorityAssignment & rcPA, int iStartLevel, vector<int> & rvectorOrder, int iWorkType, SchedCache & rcSchedCache) {
		int iTaskNum = m_pcTaskSet->getNumTasks();		
		list<int> listCandidates(rvectorOrder.begin(), rvectorOrder.end());
		for (auto cand : rvectorOrder) rcPA.unset(cand);
		int iPriority = iStartLevel;
		BitMask cLPSet(iTaskNum);
		for (int i = iTaskNum - 1; i > iPriority; i--) cLPSet.set(rcPA.getTask(i));
		for (; iPriority >= 0; iPriority--)
		{
			bool bAssigned = false;
			auto iter = listCandidates.begin();			
			for (; iter != listCandidates.end(); iter++) {
				int iCand = *iter;
				if (bAssigned = (rcSchedCache[iCand].existSubset(cLPSet).first)) {
					rcPA.setPriority(iCand, iPriority);
					cLPSet.set(iCand);
					listCandidates.erase(iter);
					break;
				}
			}

			iter = listCandidates.begin();
			for (; iter != listCandidates.end() && bAssigned == false; iter++) {
				int iCand = *iter;				
				rcPA.setPriority(iCand, iPriority);				
				auto rt = ResponseTimeAnalysis_Guan().RTTerminateOnBound(*m_pcTaskSet, iCand, rvectorRTERange[0], rcPA, rvectorRTERange[1][iCand]);
				bAssigned = rt <= rvectorRTERange[1][iCand];
				if (!bAssigned) rcPA.unset(iCand);
				else {
					listCandidates.erase(iter);
					cLPSet.set(iCand);					
					break;
				}
			}
			if (!bAssigned) break;			

			if (iWorkType == 1 && rcPA.getTask(iPriority) == iTaskIndex) {
				for (auto remaining : listCandidates) rcPA.setPriority(remaining, --iPriority);
				assert(rcPA.getSize() == iTaskNum);
				return -1;
			}
		}
		if (iPriority == -1) m_cRelaxedSchedulablePAList.push(rcPA);
		return iPriority;
	}

	void MPSchedulingUnschedCoreComputer::ComputeMaxUnschedRTELB_SchedCache(int iTaskIndex, vector<vector<ValueType>> & rvectorRTERange, pair<PriorityAssignment, int> & rcPAInfo) {
		Timer cTimer; cTimer.Start();
		numComputeMaxUnschedRTELB++;
		int iTaskNum = m_pcTaskSet->getNumTasks();
		PriorityAssignment & rcPA = rcPAInfo.first;
		ValueType cUnschedBound = rvectorRTERange[0][iTaskIndex], cSchedBound = m_pcTaskSet->getTask(iTaskIndex).C();
		SchedCache cUnschedCache, cNewUnschedCache;
		auto AddUnschedCache = [&]()->void {
			for (auto & ele : cNewUnschedCache) cUnschedCache[ele.first].insert(ele.second);			
		};
		rvectorRTERange[0][iTaskIndex] = m_pcTaskSet->getTask(iTaskIndex).C();
		if (rcPA.getPriority(iTaskIndex) != -1) {
			cTimer.Stop(); timeComputeMaxUnschedRTELB += cTimer.getElapsedTime_ms();
			return;
		}

		vector<int> vectorOrder = { iTaskIndex };
		for (int i = 0; i < iTaskNum; i++) if (rcPA.getPriority(i) <= rcPAInfo.second && i != iTaskIndex) vectorOrder.push_back(i);
		sort(vectorOrder.begin() + 1, vectorOrder.end(), [&](int lhs, int rhs)->bool {return m_pcTaskSet->getTask(lhs).D() > m_pcTaskSet->getTask(rhs).D();});
		assert(vectorOrder.size() == rcPAInfo.second + 1);

		int iStatus = ComputeMaxUnschedRTELB_isSchedulable_SchedCache(iTaskIndex, rvectorRTERange, rcPA, rcPAInfo.second, vectorOrder, 0, m_cSchedCachedForMUSConversion, cUnschedCache, cNewUnschedCache);
		if (iStatus != -1) {
			cTimer.Stop(); timeComputeMaxUnschedRTELB += cTimer.getElapsedTime_ms();
			rcPAInfo.second = iStatus;
			return;
		}
		AddUnschedCache(); cNewUnschedCache.clear();

		auto cExtractOrder = [&]()->void {
			for (int i = rcPAInfo.second; i >= 0; i--) {
				assert(rcPA.getTask(i) != -1);
				vectorOrder.at(rcPAInfo.second - i) = rcPA.getTask(i);
			}
		};
		cExtractOrder();
		AddSchedCached(rcPA, rcPA.getPriority(iTaskIndex), m_cSchedCachedForMUSConversion);

		PriorityAssignment cUnschedPA = rcPAInfo.first;
		int iUnschedFrontier = rcPAInfo.second;
		while (cUnschedBound - cSchedBound > 1) {
			rvectorRTERange[0][iTaskIndex] = (cSchedBound + cUnschedBound) >> 1;
			int iStatus = ComputeMaxUnschedRTELB_isSchedulable_SchedCache(iTaskIndex, rvectorRTERange, rcPA, rcPAInfo.second, vectorOrder, 1, m_cSchedCachedForMUSConversion, cUnschedCache, cNewUnschedCache);
			if (iStatus == -1) {
				cSchedBound = rvectorRTERange[0][iTaskIndex];
				cExtractOrder();
				AddSchedCached(rcPA, rcPA.getPriority(iTaskIndex), m_cSchedCachedForMUSConversion);
				AddUnschedCache(); cNewUnschedCache.clear();
			}
			else {
				cUnschedBound = rvectorRTERange[0][iTaskIndex];
				cUnschedPA = rcPA;
				iUnschedFrontier = iStatus;
				cNewUnschedCache.clear();
			}
		}
		assert(iUnschedFrontier != -1);
		rvectorRTERange[0][iTaskIndex] = cUnschedBound;
		for (int i = iUnschedFrontier; i >= 0; i--) cUnschedPA.unset(cUnschedPA.getTask(i));
		rcPAInfo = { cUnschedPA, iUnschedFrontier };
		cTimer.Stop();
		timeComputeMaxUnschedRTELB += cTimer.getElapsedTime_ms();
	}

	int MPSchedulingUnschedCoreComputer::ComputeMaxUnschedRTELB_isSchedulable_SchedCache(int iTaskIndex, vector<vector<ValueType>> & rvectorRTERange, PriorityAssignment & rcPA, int iStartLevel, vector<int> & rvectorOrder, int iWorkType, SchedCache & rcSchedCache, SchedCache & rcUnschedCache, SchedCache & rcNewUnschedCache) {
		int iTaskNum = m_pcTaskSet->getNumTasks();
		list<int> listCandidates(rvectorOrder.begin(), rvectorOrder.end());
		for (auto cand : rvectorOrder) rcPA.unset(cand);
		int iPriority = iStartLevel;
		BitMask cLPSet(iTaskNum);
		for (int i = iTaskNum - 1; i > iPriority; i--) cLPSet.set(rcPA.getTask(i));
		for (; iPriority >= 0; iPriority--)
		{
			bool bAssigned = false;
			auto iter = listCandidates.begin();
			for (; iter != listCandidates.end(); iter++) {
				int iCand = *iter;
				if (bAssigned = (rcSchedCache[iCand].existSubset(cLPSet).first)) {
					rcPA.setPriority(iCand, iPriority);
					cLPSet.set(iCand);
					listCandidates.erase(iter);
					break;
				}
			}

			iter = listCandidates.begin();
			for (; iter != listCandidates.end() && bAssigned == false; iter++) {
				int iCand = *iter;
				if (rcUnschedCache[iCand].existSupset(cLPSet).first) 
					continue;
				rcPA.setPriority(iCand, iPriority);
				bAssigned = DeadlineAnalysis().isSchedulable(*m_pcTaskSet, iCand, rcPA);
				if (!bAssigned) {
					auto rt = ResponseTimeAnalysis_Guan().RTTerminateOnBound(*m_pcTaskSet, iCand, rvectorRTERange[0], rcPA, rvectorRTERange[1][iCand]);
					bAssigned = rt <= rvectorRTERange[1][iCand];
				}				
				if (!bAssigned) {
					rcPA.unset(iCand);
					rcNewUnschedCache[iCand].insert(cLPSet);
				}
				else {
					listCandidates.erase(iter);
					cLPSet.set(iCand);
					break;
				}
			}
			if (!bAssigned) break;

			if (iWorkType == 1 && rcPA.getTask(iPriority) == iTaskIndex) {
				for (auto remaining : listCandidates) rcPA.setPriority(remaining, --iPriority);
				assert(rcPA.getSize() == iTaskNum);
				return -1;
			}
		}
		if (iPriority == -1) m_cRelaxedSchedulablePAList.push(rcPA);
		return iPriority;
	}

	void MPSchedulingUnschedCoreComputer::AddSchedCached(PriorityAssignment & rcPA, int iLB, SchedCache & rcSchedCachedForMUSConversion){
		BitMask rcLPSet(m_pcTaskSet->getNumTasks());
		for (int i = m_pcTaskSet->getNumTasks() - 1; i > iLB; i--) rcLPSet.set(rcPA.getTask(i));
		for (int i = iLB; i >= 0; i--) {						
			rcSchedCachedForMUSConversion[rcPA.getTask(i)].insert(rcLPSet);
			rcLPSet.set(rcPA.getTask(i));
		}
	}

	//ExhaustiveSearch_GuanRTA
	bool ExhaustiveSearch_GuanRTA::recur(int iLevel, MultiProcessorSystem & rcSystem, PriorityAssignment & rcPA) {
		if (iLevel == rcSystem.getNumTasks()) {
			vector<ValueType> vecRT(rcSystem.getNumTasks(), 0);
			if (count++ % 1000 == 0) cout << "\rExamined: " << (count) / 1000 << "k    ";
			if (ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(rcSystem, rcPA, vecRT)) return true;
			return false;
		}
		for (int i = 0; i < rcSystem.getNumTasks(); i++) {
			if (rcPA.getPriority(i) != -1) continue;
			rcPA.setPriority(i, iLevel);
			if (recur(iLevel + 1, rcSystem, rcPA)) return true;
			rcPA.unset(i);
		}
		return false;
	}

	pair<bool, PriorityAssignment> ExhaustiveSearch_GuanRTA::operator () (MultiProcessorSystem & rcSystem) {
		PriorityAssignment cPA(rcSystem);
		count = 0;
		bool status = recur(0, rcSystem, cPA);
		cout << endl;
		return { status, cPA };
	}
	

#ifdef MILP		
	MPSchedulingFR::Configuration::Configuration() {
		m_iMaximizeMinimumSlackObjective = 0;
		m_iUseAsyncTest = 1;
	}

	MPSchedulingFR::MPSchedulingFR(SystemType & rcSystem)
		:GeneralFocusRefinement_MonoInt(rcSystem) {

	}

	void MPSchedulingFR::CreateSchedCondTypeVars(const IloEnv & rcEnv, SchedCondTypeVars & rcVars) {
		int iTaskNum = m_pcTaskSet->getNumTasks();
		char axBuffer[512] = { 0 };
		for (int i = 0; i < iTaskNum; i++) {
			sprintf(axBuffer, "DUB(%d)", i);
			rcVars[SchedCondType(i, 0)] = IloNumVar(rcEnv, m_pcTaskSet->getTask(i).C(), m_pcTaskSet->getTask(i).D(), axBuffer);
			sprintf(axBuffer, "DLB(%d)", i);
			rcVars[SchedCondType(i, 1)] = IloNumVar(rcEnv, m_pcTaskSet->getTask(i).C(), m_pcTaskSet->getTask(i).D(), axBuffer);
		}
	}

	void MPSchedulingFR::GenerateLBUBEqualConst(const IloEnv & rcEnv, SchedCondTypeVars & rcVars, IloRangeArray & rcConst) {
		int iTaskNum = m_pcTaskSet->getNumTasks();
		for (int i = 0; i < iTaskNum; i++) {
			assert(rcVars.count(SchedCondType(i, 0)));
			assert(rcVars.count(SchedCondType(i, 1)));
			rcConst.add(rcVars[SchedCondType(i, 0)] - rcVars[SchedCondType(i, 1)] == 0);
		}
	}

	int MPSchedulingFR::SolveFR(int iDisplay, SchedConds & rcSchedCondFixed, UnschedCores & rcGenUC, double dObjLB, double dObjUB, double dTimeout) {
		try
		{
			IloEnv cSolveEnv;
			IloRangeArray cConst(cSolveEnv);
			SchedCondTypeVars cSchedCondTypeVars;
			SchedCondVars cSchedCondVars;
			IloNumVar cDummyInteger(cSolveEnv, 0.0, 1.0, IloNumVar::Int);
			CreateSchedCondTypeVars(cSolveEnv, cSchedCondTypeVars);
			CreateSchedCondVars(cSolveEnv, rcGenUC, cSchedCondVars);
			GenSchedCondTypeConst(cSolveEnv, cSchedCondTypeVars, cSchedCondVars, cConst);
			GenUnschedCoreConst(cSolveEnv, rcGenUC, cSchedCondVars, cConst);
			GenerateLBUBEqualConst(cSolveEnv, cSchedCondTypeVars, cConst);						
			IloObjective cObjective = (m_cConfiguration.m_iMaximizeMinimumSlackObjective == 1) ?
				getMaximizeMinimumSlackObjective(cSolveEnv, cSchedCondTypeVars, cConst) : IloMinimize(cSolveEnv, 0);
			IloModel cModel(cSolveEnv);
			cModel.add(cConst);			
			cModel.add(cDummyInteger);		
			cModel.add(cObjective);			
			IloCplex cSolver(cModel);
			EnableSolverDisp(iDisplay, cSolver);
			setSolverParam(cSolver);

			cSolver.setParam(IloCplex::Param::TimeLimit, dTimeout);
			int iRetStatus = 0;

			m_cSchedCondsConfig.clear();
			m_dequeSolutionSchedConds.clear();
			m_iBestSolFeasible = 0;
			Timer cSolveTimer;cSolveTimer.Start();
			m_dSolveWallTime = cSolver.getCplexTime();
			bool bStatus = cSolver.solve();			
			cSolveTimer.Stop();
			m_dSolveCPUTime = cSolveTimer.getElapsedTime_ms();
			m_dSolveWallTime = cSolver.getCplexTime() - m_dSolveWallTime;
			m_cCplexStatus = cSolver.getCplexStatus();

			m_dequeSolutionSchedConds.clear();
			m_cSchedCondsConfig.clear();			
			if (bStatus) {
				ExtractAllSolution(cSolver, cSchedCondTypeVars);
				MarshalAllSolution(m_dequeSolutionSchedConds);
			}

			if ((cSolver.getCplexStatus() == cSolver.Optimal) || (cSolver.getCplexStatus() == cSolver.OptimalTol))
				iRetStatus = 1;
			else if (cSolver.getCplexStatus() == cSolver.AbortTimeLim)
				iRetStatus = -1;
			else if (cSolver.getCplexStatus() == cSolver.AbortUser)
				iRetStatus = -2;
			else if (cSolver.getCplexStatus() == cSolver.Infeasible)
				iRetStatus = 0;
			else
				iRetStatus = 0;
			cSolveEnv.end();
			return iRetStatus;

		}
		catch (IloException& e)
		{
			cerr << "Concert exception caught: " << e << endl;
			while (1);
		}
	}

	bool MPSchedulingFR::AsynchronousTestControl::empty() {
		m_cLock.lock();
		bool bStatus = m_dequeFound.empty();
		m_cLock.unlock();
		return bStatus;
	}

	void MPSchedulingFR::AsynchronousTestControl::push(PriorityAssignment & rcPA) {
		m_cLock.lock();	
		m_dequeFound.push_back(rcPA);
		m_cLock.unlock();		
	}

	void MPSchedulingFR::AsynchronousTestControl::terminate() {
		m_cLock.lock();		
		m_bTerminate = true;
		m_cLock.unlock();
	}

	bool MPSchedulingFR::AsynchronousTestControl::isTerminate() {
		m_cLock.lock();
		bool bRet = m_bTerminate;
		m_cLock.unlock();
		return bRet;
	}

	PriorityAssignment MPSchedulingFR::AsynchronousTestControl::top() {
		m_cLock.lock();
		assert(m_dequeFound.empty() == false);
		auto cRet = m_dequeFound.front();
		m_cLock.unlock();
		return cRet;
	}
	
	void MPSchedulingFR::StartAsynchronousTest(int iThreadNum) {
		m_vectorAsynchronousTestThread = vector<thread>(iThreadNum);
		for (int i = 0; i < iThreadNum; i++) {
			m_vectorAsynchronousTestThread[i] = thread([&]()->void {
				while (m_cAsyncTestControl.isTerminate() == false && m_cAsyncTestControl.empty()) {
					auto cSample = m_cUnschedCoreComputer.getRelaxedSchedulablePAList().pop();
					if (cSample.empty()) {
						Sleep(1000);
						continue;
					}
					auto cResult = ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(*m_pcTaskSet, cSample.front());
					if (cResult.first) m_cAsyncTestControl.push(cSample.front());
				}
			});
		}
	}

	void MPSchedulingFR::CloseAsynchronousTest() {
		m_cAsyncTestControl.terminate();
		for (auto & ele : m_vectorAsynchronousTestThread) ele.join();
	}

	int MPSchedulingFR::FocusRefinement(int iDisplay, double dTimeout)
	{
		assert(m_pcTaskSet);		
		if(m_cConfiguration.m_iUseAsyncTest) StartAsynchronousTest(1);
		m_dStartSeconds = getWallTimeSecond();
		int iRetStatus = 0;							
		Timer cTimer; cTimer.Start();		
		int iIterationCount = 0;
		double dUCComputeTime = 0;
		bool bCanExit = false;		
		while (!bCanExit){
			bCanExit = true;			
			if (iDisplay == 2)	cout << "ILP query...." << endl;
			PriorityAssignment cPA(*m_pcTaskSet);
			int iStatus = 0;
			if (m_cAsyncTestControl.empty() == false) {
				cPA = m_cAsyncTestControl.top();
				auto cResult = ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(*m_pcTaskSet, cPA);
				for (int i = 0; i < m_pcTaskSet->getNumTasks(); i++) {
					m_cSchedCondsConfig[SchedCondType(i, 0)] = cResult.second[i];
					m_cSchedCondsConfig[SchedCondType(i, 1)] = cResult.second[i];
				}
				iStatus = 1;
			}
			else
				iStatus = SolveFR(m_iSubILPDisplay, m_cSchedCondsFixed, m_cUnschedCores, m_dObjLB, 1e74, 1e74);

			m_dequeIterationLogs.push_back(IterationStatus());
			IterationStatus & tagStatus = m_dequeIterationLogs.back();
			tagStatus.dTime = (double)(getWallTimeSecond() - m_dStartSeconds); ;
			tagStatus.cSchedCondsConfig = m_cSchedCondsConfig;
			tagStatus.iIterationN = iIterationCount++;
			tagStatus.dObjective = m_dObjective;
			tagStatus.dBestFeasible = m_dObjUB;
			tagStatus.dObjLB = m_dObjLB;
			tagStatus.iTotalCores = m_cUnschedCores.size();			
			tagStatus.dILPTime = m_dSolveWallTime;
			cTimer.Stop();
			tagStatus.dCPUTime = cTimer.getElapsedTime_ms();

			if (iStatus == 0) {
				tagStatus.enumState = IterationState::Infeasible;
				iRetStatus = 0;
			}
			else if (tagStatus.dTime > dTimeout) {
				tagStatus.enumState = Timeout;
				iRetStatus = -1;
			}
			else {
				bool bSchedulable = m_cUnschedCoreComputer.IsSchedulable(m_cSchedCondsConfig, m_cSchedCondsFixed, cPA);
				if (!bSchedulable) {
					UnschedCores cThisCores;
					double dUCSeconds = getWallTimeSecond();
					if (iDisplay == 2) cout << "Unched-core computation..." << endl;
					m_cUnschedCoreComputer.ComputeUnschedCores(m_cSchedCondsConfig, m_cSchedCondsFixed, cThisCores, m_iCorePerIter);
					for (auto ele : cThisCores) m_cUnschedCores.push_back(ele);
					tagStatus.dUCComputeTime = getWallTimeSecond() - dUCSeconds;
					tagStatus.enumState = IterationState::Search;
					tagStatus.cThisUnschedCores = cThisCores;
					bCanExit = false;
				}
				else {
					m_cPriorityAssignmentSol = cPA;
					tagStatus.enumState = Optimal;
					if (ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(*m_pcTaskSet, cPA).first == false) {
						cout << "Response time verification failed? " << endl;
					}
					assert(ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(*m_pcTaskSet, cPA).first);
					iRetStatus = 1;
				}
			}																											
			//Display Status
			if (iDisplay) PrintIterationStatusSimple(m_dequeIterationLogs.back());
			PrintIterStatusToLog(m_dequeIterationLogs.size() - 1);			
		}
		if (iDisplay) cout << endl;
		ComputeStatistic();		
		if (iDisplay == 2) cout << "Closing asynchronous test" << endl;
		if (m_cConfiguration.m_iUseAsyncTest) CloseAsynchronousTest();
		if (iDisplay == 2) cout << "Asynchronous test closed" << endl;
		return iRetStatus;
	}


	IloObjective MPSchedulingFR::getMaximizeMinimumSlackObjective(const IloEnv & rcEnv, SchedCondTypeVars & rcVars, IloRangeArray & rcConst) {
		int iTaskNum = m_pcTaskSet->getNumTasks();
		IloNumVar cMinSlack(rcEnv, "MinSlack");
		for (int i = 0; i < iTaskNum; i++) {
			rcConst.add(cMinSlack - (m_pcTaskSet->getTask(i).D() - rcVars[SchedCondType(i, 0)]) <= 0);
		}		
		return IloMinimize(rcEnv, -cMinSlack);
	}
	
#endif
	extern double dRTTerminateOnBoundTime;
	void TestOptimizedUnschedCoreComputer() {
		MultiProcessorSystem cSystem;
		//cSystem.Read("F:\\C++Project\\MPSchedulingPA\\MPSchedulingPA\\MPSchedulingPA\\ExpSpace\\TestMPSystem.txt");
		cSystem.Read("F:\\C++Project\\MPSchedulingPA\\MPSchedulingPA\\MPSchedulingPA\\ExpSpace\\System.txt");
		string stringSchedConds = "({0, 0} : 10534), ({0, 1} : 10534), ({1, 1} : 48), ({2, 0} : 26074), ({2, 1} : 26074), ({3, 0} : 511), ({3, 1} : 511), ({4, 0} : 53754), ({4, 1} : 53754), ({5, 0} : 941), ({5, 1} : 941), ({6, 0} : 3341), ({6, 1} : 3341), ({7, 0} : 2332), ({7, 1} : 2332), ({8, 0} : 36859), ({9, 0} : 5387), ({9, 1} : 5387), ({10, 0} : 79943), ({10, 1} : 79943), ({11, 0} : 1533), ({11, 1} : 1533), ({12, 0} : 195), ({13, 0} : 2604), ({13, 1} : 2604), ({14, 0} : 4224), ({14, 1} : 4224), ({15, 0} : 318), ({15, 1} : 318), ({16, 0} : 5886), ({16, 1} : 5886), ({17, 0} : 501), ({17, 1} : 501), ({18, 0} : 2158), ({18, 1} : 2158), ({19, 0} : 14408), ({19, 1} : 14408), ({20, 0} : 15417), ({20, 1} : 15417), ({21, 0} : 68), ({21, 1} : 68), ({22, 0} : 48863), ({22, 1} : 48863), ({23, 0} : 59664), ({23, 1} : 59664), ({24, 0} : 17374), ({24, 1} : 17374), ({25, 0} : 63), ({25, 1} : 63), ({26, 0} : 1337), ({26, 1} : 1337), ({27, 0} : 78610), ({27, 1} : 78610), ({28, 0} : 432), ({29, 0} : 59424), ({29, 1} : 59424), ({30, 0} : 322), ({30, 1} : 322), ({31, 0} : 75), ({31, 1} : 75), ({32, 0} : 2330), ({32, 1} : 2330), ({33, 0} : 14601), ({33, 1} : 14601), ({34, 0} : 1380), ({34, 1} : 1380), ({35, 0} : 3024), ({35, 1} : 3024), ({36, 0} : 98808), ({37, 0} : 6222), ({37, 1} : 6222), ({38, 0} : 199), ({38, 1} : 199), ({39, 0} : 93), ({39, 1} : 93), ({40, 0} : 5611), ({40, 1} : 5611), ({41, 0} : 101), ({42, 0} : 1137), ({42, 1} : 1137), ({43, 0} : 1454), ({43, 1} : 1454), ({44, 0} : 17071), ({44, 1} : 17071), ({45, 0} : 65429), ({46, 0} : 14483), ({46, 1} : 14483), ({47, 0} : 5158), ({47, 1} : 5158), ({48, 0} : 16045), ({48, 1} : 16045), ({49, 0} : 1267), ({49, 1} : 1267)";//Iteration 77
		stringSchedConds = "({0, 1} : 3100), ({1, 1} : 48), ({2, 1} : 2329), ({3, 1} : 74), ({4, 1} : 1664), ({5, 1} : 77), ({6, 1} : 1169), ({7, 1} : 272), ({8, 1} : 199), ({9, 1} : 415), ({10, 1} : 6844), ({11, 1} : 265), ({12, 1} : 1), ({13, 1} : 403), ({14, 1} : 162), ({15, 1} : 11), ({16, 1} : 905), ({17, 1} : 145), ({18, 1} : 182), ({19, 1} : 482), ({20, 1} : 3973), ({21, 1} : 12), ({22, 1} : 2239), ({23, 1} : 1763), ({24, 1} : 2228), ({25, 1} : 5), ({26, 1} : 627), ({27, 1} : 3063), ({28, 1} : 1), ({29, 1} : 6558), ({30, 1} : 8), ({31, 1} : 21), ({32, 1} : 112), ({33, 1} : 867), ({34, 1} : 268), ({35, 1} : 58), ({36, 1} : 99), ({37, 1} : 218), ({38, 1} : 4), ({39, 1} : 13), ({40, 1} : 680), ({41, 1} : 1), ({42, 1} : 49), ({43, 1} : 179), ({44, 1} : 1795), ({45, 1} : 1068), ({46, 1} : 722), ({47, 1} : 119), ({48, 1} : 1518), ({49, 1} : 19)";
		stringSchedConds = "({0, 1} : 778), ({1, 1} : 2865), ({2, 1} : 8), ({3, 1} : 586), ({4, 1} : 57), ({5, 1} : 41), ({6, 1} : 3873), ({7, 1} : 27), ({8, 1} : 39), ({9, 1} : 64), ({10, 1} : 156), ({11, 1} : 591), ({12, 1} : 141), ({13, 1} : 1123), ({14, 1} : 2391), ({15, 1} : 58), ({16, 1} : 1), ({17, 1} : 55165), ({18, 1} : 4417), ({19, 1} : 131)";
		stringSchedConds = "({0, 1} : 520), ({1, 1} : 1616), ({2, 1} : 1198), ({3, 1} : 276), ({4, 1} : 453), ({5, 1} : 277), ({6, 1} : 8952), ({7, 1} : 192), ({8, 1} : 12), ({9, 1} : 83), ({10, 1} : 8094), ({11, 1} : 60), ({12, 1} : 26), ({13, 1} : 53), ({14, 1} : 174), ({15, 1} : 68), ({16, 1} : 30), ({17, 1} : 5), ({18, 1} : 45), ({19, 1} : 2982), ({20, 1} : 14), ({21, 1} : 23), ({22, 1} : 77), ({23, 1} : 1402), ({24, 1} : 2666), ({25, 1} : 2), ({26, 1} : 275), ({27, 1} : 4207), ({28, 1} : 15498), ({29, 1} : 152), ({30, 1} : 4251), ({31, 1} : 4402), ({32, 1} : 26), ({33, 1} : 98), ({34, 1} : 27817), ({35, 1} : 701), ({36, 0} : 7222), ({37, 1} : 1376), ({38, 0} : 57433), ({39, 0} : 6599), ({40, 0} : 75091), ({41, 1} : 1), ({42, 1} : 40), ({43, 0} : 26380), ({44, 0} : 6243), ({45, 1} : 23), ({46, 0} : 1379), ({47, 1} : 189), ({48, 0} : 86914), ({49, 0} : 14027), ({50, 0} : 6928), ({51, 0} : 133), ({52, 0} : 858), ({53, 0} : 138), ({54, 0} : 2185), ({55, 0} : 670), ({56, 0} : 164), ({57, 0} : 318), ({58, 0} : 338), ({59, 0} : 113), ({60, 1} : 4607), ({61, 1} : 226), ({62, 0} : 35522), ({63, 0} : 12844), ({64, 0} : 4338), ({65, 1} : 47), ({66, 1} : 22622), ({67, 0} : 15602), ({68, 1} : 18), ({69, 1} : 3), ({70, 1} : 820), ({71, 0} : 30344), ({72, 1} : 178), ({73, 0} : 44557), ({74, 0} : 2681), ({75, 1} : 2789), ({76, 0} : 96024), ({77, 0} : 78267), ({78, 0} : 40429), ({79, 0} : 40106)";
		//stringSchedConds = "({0, 1} : 520), ({1, 1} : 1616), ({2, 1} : 1198), ({3, 1} : 276), ({4, 1} : 453), ({5, 1} : 277), ({6, 1} : 8952), ({7, 1} : 192), ({8, 1} : 12), ({9, 1} : 83), ({10, 1} : 8094), ({11, 1} : 60), ({12, 1} : 26), ({13, 1} : 53), ({14, 1} : 174), ({15, 1} : 68), ({16, 1} : 30), ({17, 1} : 5), ({18, 1} : 45), ({19, 1} : 2982), ({20, 1} : 14), ({21, 1} : 23), ({22, 1} : 77), ({23, 1} : 1402), ({24, 1} : 2666), ({25, 1} : 2), ({26, 1} : 275), ({27, 1} : 4207), ({28, 1} : 15498), ({29, 1} : 152), ({30, 1} : 4251), ({31, 1} : 4402), ({32, 1} : 26), ({33, 1} : 98), ({34, 1} : 27817), ({35, 1} : 701), ({36, 0} : 7222), ({37, 1} : 1376), ({38, 0} : 57433), ({39, 1} : 2356), ({40, 0} : 75091), ({41, 1} : 1), ({42, 1} : 40), ({43, 1} : 3306), ({44, 1} : 616), ({45, 0} : 448), ({46, 0} : 219), ({46, 1} : 219), ({47, 0} : 583), ({47, 1} : 583), ({48, 0} : 86914), ({49, 1} : 8680), ({50, 0} : 1712), ({50, 1} : 1712), ({51, 0} : 75), ({51, 1} : 75), ({52, 0} : 278), ({52, 1} : 278), ({53, 0} : 86), ({53, 1} : 86), ({54, 1} : 424), ({55, 0} : 166), ({55, 1} : 166), ({56, 1} : 67), ({57, 0} : 318), ({58, 1} : 26), ({59, 0} : 33), ({59, 1} : 33), ({60, 1} : 4607), ({61, 0} : 336), ({61, 1} : 336), ({62, 0} : 35522), ({63, 0} : 12844), ({64, 0} : 4338), ({65, 1} : 47), ({66, 1} : 22622), ({67, 0} : 15602), ({68, 1} : 18), ({69, 1} : 3), ({70, 1} : 820), ({71, 0} : 30344), ({72, 1} : 178), ({73, 0} : 44557), ({74, 0} : 1775), ({74, 1} : 1775), ({75, 1} : 2789), ({76, 0} : 96024), ({77, 0} : 78267), ({78, 0} : 40429), ({79, 0} : 40106)";
		MPSchedulingUnschedCoreComputer cUCComp(cSystem);
		MPSchedulingUnschedCoreComputer::SchedConds cFlex, cDummy;
		Timer cTimer; cTimer.Start();
		cUCComp.ReadSchedCondsFromStr(cFlex, stringSchedConds.data());
		MPSchedulingUnschedCoreComputer::UnschedCores cCores;
		cUCComp.ComputeUnschedCores(cFlex, cDummy, cCores, 5);
		cTimer.Stop();
		for (auto ele : cCores) {
			cUCComp.PrintSchedConds(ele, cout); cout << endl;
		}
		cout << cTimer.getElapsedTime_ms() << endl;
		cout << dRTTerminateOnBoundTime << endl;
		cout << timeComputeMaxUnschedRTELB << endl;
		cout << numComputeMaxUnschedRTELB << endl;
	}
#ifdef MILP
	void DebugFR() {
		MultiProcessorSystem cSystem;
		//cSystem.Read("F:\\C++Project\\MPSchedulingPA\\MPSchedulingPA\\MPSchedulingPA\\ExpSpace\\M4N20\\Util0.500000\\Case0\\System.txt");
		//cSystem.Read("F:\\C++Project\\MPSchedulingPA\\MPSchedulingPA\\MPSchedulingPA\\ExpSpace\\M4N20\\Util4.000000\\Case2\\System.txt");
		//cSystem.Read("F:\\C++Project\\MPSchedulingPA\\MPSchedulingPA\\MPSchedulingPA\\ExpSpace\\M8N40\\System6.20-450.txt");
		cSystem.Read("System.txt");
		//string stringSchedConds = "({0, 1} : 778), ({1, 1} : 2865), ({2, 1} : 8), ({3, 1} : 586), ({4, 1} : 57), ({5, 1} : 41), ({6, 1} : 3873), ({7, 1} : 27), ({8, 1} : 39), ({9, 1} : 64), ({10, 1} : 156), ({11, 1} : 591), ({12, 1} : 141), ({13, 1} : 1123), ({14, 1} : 2391), ({15, 1} : 58), ({16, 1} : 1), ({17, 1} : 55165), ({18, 1} : 4417), ({19, 1} : 131)";
		MPSchedulingUnschedCoreComputer cUnschedCoreComputer(cSystem);
		
		//DA
		PriorityAssignment cDAPA(cSystem);
		int iDAResult = DeadlineAnalysis().Audsley(cSystem, cDAPA, 0);
		cout << "DA Result: " << iDAResult << endl;
		auto cResult = ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(cSystem, cDAPA);

		MPSchedulingFR cFR(cSystem);
		cFR.setCorePerIter(5);
		//cFR.setSubILPDisplay(2);
		//cFR.getConfiguration().m_iMaximizeMinimumSlackObjective = 1;
		cFR.setLogFile("DebugFRLog.txt");
		int status = cFR.FocusRefinement(1);
		cout << "Debug end" << endl;

	}
#endif

	void TestMPSchedulingFR() {
		RandomSystemGenerator cGenerator;
		MultiProcessorSystem cSystem;
#if 0
		double dUtil = 1.0;
		cin >> dUtil;
		cSystem = cGenerator.GenerateSystem(50, 4, dUtil, 100, 100000);
		cSystem.Write("TestMPSystem.txt");
#else
		cSystem.Read("TestMPSystem.txt");
#endif
		MPSchedulingUnschedCoreComputer::SchedConds cDummy;
		MPSchedulingUnschedCoreComputer cUnschedCoreComputer(cSystem);
		cout << cUnschedCoreComputer.IsSchedulable(cDummy, cDummy) << endl;
		cout << "System Utilization: " << cSystem.TotalUtil() << endl;

		//DkC Priority Assignment
		cout << "DkC PA" << endl;
		PriorityAssignment cDkCPA; cDkCPA.GenerateDkCPA(cSystem);
		auto DkCResult = ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(cSystem, cDkCPA);
		cout << "Status: " << DkCResult.first << endl;
		for (int i = 0; i < cSystem.getNumTasks(); i++) {
			//cout << "Task " << i << ": " << DkCResult.second[i] << ' ' << cSystem.getTask(i).D() - DkCResult.second[i] << endl;
		}

		//DM Priority Assignment
		cout << "DM PA" << endl;
		PriorityAssignment cDMPA; cDMPA.GenerateDMPA(cSystem);
		auto DMResult = ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(cSystem, cDMPA);
		cout << "Status: " << DMResult.first << endl;
		for (int i = 0; i < cSystem.getNumTasks(); i++) {
			//cout << "Task " << i << ": " << DMResult.second[i] << ' ' << cSystem.getTask(i).D() - DMResult.second[i] << endl;
		}

#ifdef MILP
		MPSchedulingFR cFR(cSystem);
		cFR.setLogFile("TestSystemFRLog.txt");
		//cFR.setSubILPDisplay(2);
		cFR.setCorePerIter(5);
		int status = cFR.FocusRefinement(1);	
		cFR.GenerateStatisticFile("TestSystemStatistic.txt");
		if (status == 1) {
			auto cSchedConds = cFR.getOptimalSchedCondConfig();
			vector<ValueType> vecRTE(cSystem.getNumTasks());
			for (int i = 0; i < cSystem.getNumTasks(); i++) {
				auto cType = MPSchedulingUnschedCoreComputer::SchedCondType(i, 0);
				vecRTE[i] = cSchedConds.count(cType) ? cSchedConds[cType] : cSystem.getTask(i).C();
			}

			vector<ValueType> vecGuanRTE(cSystem.getNumTasks(), 0.0);
			auto cPA = cFR.getSolPA();
			cPA.Write("TestMPSystemPA.txt");
			ResponseTimeAnalysis_Guan().ComputeAllRTGivenPA(cSystem, cPA, vecGuanRTE);
			for (int i = 0; i < cSystem.getNumTasks(); i++) {
				cout << "Task " << i << ": " << vecRTE[i] << " " << vecGuanRTE[i] << ' ' << cSystem.getTask(i).D() - vecGuanRTE[i] << endl;
			}						
		}			
#endif		
		cout << "isSchedulableTime: " << isSchedulableTime << endl;
		cout << "timeComputeMaxUnschedRTELB: " << timeComputeMaxUnschedRTELB << endl;
	}
}