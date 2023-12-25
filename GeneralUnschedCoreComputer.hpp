#pragma once
#include <stddef.h>
#include <stdint.h>
#include <iostream>
using namespace std;
#include <vector>
#include <map>
#include <time.h>
#include <deque>
#include <sstream>
#include "MultiProcessorSystem.h"
#include <assert.h>
#include "Util.h"

using namespace MPScheduling;
template<class System, class SchedCondKey, class SchedCondData = int64_t>
class GeneralUnschedCoreComputer_MonoInt
{
public:

	typedef SchedCondKey SchedCondType;
	typedef SchedCondData SchedCondContent;
	typedef std::pair<SchedCondType, SchedCondContent> SchedCond;
	typedef map<SchedCondType, SchedCondContent> SchedConds;
	typedef SchedConds UnschedCore;
	typedef deque<SchedConds> SchedCondsArray;
	typedef SchedCondsArray UnschedCores;
	typedef System SystemType;

	enum Monotonicity
	{
		Up,
		Down		
	};

protected:
	SystemType * m_pcTaskSet;

	//Enforced Partial Order 
	typedef vector<int> IntArray1D;
	typedef vector<IntArray1D> IntArray2D;
	IntArray2D m_cEnforcedPartialOrder;
public:
	GeneralUnschedCoreComputer_MonoInt()
	{

	}

	GeneralUnschedCoreComputer_MonoInt(SystemType & rcSystem)
	{
		m_pcTaskSet = &rcSystem;
		InitiailizePartialOrderTable();
	}

	~GeneralUnschedCoreComputer_MonoInt()
	{

	}

	virtual bool IsSchedulable()
	{
		SchedConds cDummy;
		return IsSchedulable(cDummy, cDummy);
	}

	virtual bool IsSchedulable(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed) = 0;

	virtual bool IsSchedulable(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed, PriorityAssignment & rcPriorityAssignment) = 0;

	int ComputeUnschedCores(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedConstFixed, SchedCondsArray & rcUnschedCores, int iLimit, double dTimeout = 1e74)
	{
		assert(m_pcTaskSet);
		PreComputeUnschedCores(rcSchedCondsFlex, rcSchedConstFixed, rcUnschedCores, iLimit, dTimeout);
		int iIniSize = rcUnschedCores.size();
		int iNextLevel = 0;
		int iStatus = Recur(iNextLevel, iLimit, rcSchedCondsFlex, rcSchedConstFixed, rcUnschedCores, time(NULL), dTimeout);
		while (iStatus == -1)
		{
			iNextLevel++;
			iStatus = Recur(iNextLevel, iLimit, rcSchedCondsFlex, rcSchedConstFixed, rcUnschedCores, time(NULL), dTimeout);
		}
		return rcUnschedCores.size() - iIniSize;

	}

	void PrintSchedCond(const SchedCond & rcSchedCond, ostream & os)
	{
		os << "(";
		os << (SchedCondType &)rcSchedCond.first;
		os << " : ";
		os << rcSchedCond.second;
		os << ")";
	}

	void PrintSchedConds(SchedConds & rcSchedConds, ostream & os)
	{
		for (typename SchedConds::iterator iter = rcSchedConds.begin(); iter != rcSchedConds.end(); iter++)
		{
			if (iter != rcSchedConds.begin())
			{
				os << ", ";
			}
			PrintSchedCond(*iter, os);
		}
	}

	virtual SchedCondContent SchedCondUpperBound(const SchedCondType & rcSchedCondType)
	{
		cout << "In SchedCondUpperBound. Base class don't assume any specific sched cond. Don't use" << endl;
		while (1);
		SchedCondContent cDummy;
		return cDummy;
	}

	virtual SchedCondContent SchedCondLowerBound(const SchedCondType & rcSchedCondType)
	{
		cout << "In SchedCondLowerBound. Base class don't assume any specific sched cond. Don't use" << endl;
		while (1);
		SchedCondContent cDummy;
		return cDummy;
	}

	virtual Monotonicity SchedCondMonotonicity(const SchedCondType & rcSchedCondType)
	{
		cout << "In SchedCondLowerBound. Base class don't assume any specific sched cond. Don't use" << endl;
		while (1);
		return Up;
	}

	bool IsSchedCondRedundant(const SchedCond & rcSchedCond)
	{
		Monotonicity cMonotonicity = SchedCondMonotonicity(rcSchedCond.first);
		if (cMonotonicity == Up)
		{
			return rcSchedCond.second == SchedCondUpperBound(rcSchedCond.first);
		}
		else if (cMonotonicity == Down)
		{
			return rcSchedCond.second == SchedCondLowerBound(rcSchedCond.first);
		}
		else
		{
			cout << "Unexpected Monotonicity Type" << endl;
			while (1);
		}
		return false;
	}

	int CompareSchedCond(const SchedCond & rcLHS, const SchedCond & rcRHS)
	{
		//0: equal; -1: stricter; 1: more relaxed
		if ((rcLHS.first < rcRHS.first) || (rcRHS.first < rcLHS.first))
		{
			//incomparable
			return -2;
		}

		Monotonicity cMonotonicity = SchedCondMonotonicity(rcLHS.first);
		if (cMonotonicity == Up)
		{
			if (DOUBLE_EQUAL(rcLHS.second, rcRHS.second, 1e-7))
			{
				return 0;
			}
			else if (rcLHS.second < rcRHS.second)
			{
				return -1;
			}
			else
			{
				return 1;
			}
			
		}
		else if (cMonotonicity == Down)
		{
			if (DOUBLE_EQUAL(rcLHS.second, rcRHS.second, 1e-7))
			{
				return 0;
			}
			else if (rcLHS.second < rcRHS.second)
			{
				return 1;
			}
			else
			{
				return -1;
			}
		}
		else
		{
			cout << "Unexpected Monotonicity Type" << endl;
			while (1);
		}
		return -2;
	}

	int CompareSchedConds(SchedConds & rcLHS, SchedConds & rcRHS)
	{
		typename SchedConds::iterator iterLHS = rcLHS.begin();
		typename SchedConds::iterator iterRHS = rcRHS.begin();
		int iStricter = 1;
		int iRelaxed = 1;
		for (; iterLHS != rcLHS.end() && iterRHS != rcRHS.end();)
		{
			if (iterLHS->first < iterRHS->first)
			{				
				//LHS has potentially more restriction.
				if ((iRelaxed == 1) && (!IsSchedCondRedundant(*iterLHS)))
				{
					//The extra restriction is not redundant
					iRelaxed = 0;//can only be stricter					
				}
				iterLHS++;
			}
			else if (iterRHS->first < iterLHS->first)
			{
				//RHS has potentially more restriction.
				if ((iStricter == 1) && (!IsSchedCondRedundant(*iterRHS)))
				{
					//The extra restriction is not redundant
					iStricter = 0;//can only be more relaxed
				}
				iterRHS++;
			}
			else
			{
				int iCmpStatus = CompareSchedCond(*iterLHS, *iterRHS);
				if (iCmpStatus == 0)
				{

				}
				else if (iCmpStatus == -1)
				{
					//LHS stricter					
					iRelaxed = 0;
				}
				else if (iCmpStatus == 1)
				{
					//LHS relaxed					
					iStricter = 0;
				}
				iterLHS++;
				iterRHS++;
			}

			if (iStricter == 0 && iRelaxed == 0)
			{
				//Incomparable
				return -2;
			}
		}

		if (iterLHS == rcLHS.end())
		{
			//RHS has more conditions.
			if (iStricter == 1)
			{
				for (; iterRHS != rcRHS.end(); iterRHS++)
				{
					if (!IsSchedCondRedundant(*iterRHS))
					{
						iStricter = 0;
						break;
					}
				}
			}			
		}

		if (iterRHS == rcRHS.end())
		{
			//LHS has more conditions.
			if (iRelaxed == 1)
			{
				for (; iterLHS != rcLHS.end(); iterLHS++)
				{
					if (!IsSchedCondRedundant(*iterLHS))
					{
						iRelaxed = 0;
						break;
					}
				}
			}
			
		}

		if (iRelaxed == 1 && iStricter == 1)
		{
			return 0;
		}
		else if (iRelaxed == 1)
		{
			return 1;
		}
		else if (iStricter == 1)
		{
			return -1;
		}
		else
		{
			return -2;
		}
	}	

	void ReadSchedCondsFromStr(SchedConds & rcSchedConds, const char axBuffer[])
	{
		stringstream ss(axBuffer);
		while (ss.eof() == false)
		{
			SchedCond cNewSchedCond;
			char xChar = 0;
			ss >> xChar;
			if (xChar == ',') continue;
			if (xChar == 0)	break;
			if (xChar != '(')
			{
				cout << "Unexpected Format. Do not proceed." << endl;
				while (1);
			}
			ss >> cNewSchedCond.first;
			ss >> xChar;
			if (xChar != ':')
			{
				cout << "Unexpected Format. Do not proceed." << endl;
				while (1);
			}
			ss >> cNewSchedCond.second;
			ss >> xChar;
			if (xChar != ')')
			{
				cout << "Unexpected Format. Do not proceed." << endl;
				while (1);
			}
			rcSchedConds.insert(cNewSchedCond);
		}
	}

	SchedConds ConjunctSchedConds(SchedConds & rcLHS, SchedConds & rcRHS)
	{
		SchedConds cRet;
		typename SchedConds::iterator iterLHS = rcLHS.begin();
		typename SchedConds::iterator iterRHS = rcRHS.begin();
		int iStricter = 1;
		int iRelaxed = 1;
		for (; iterLHS != rcLHS.end() && iterRHS != rcRHS.end();)
		{
			if (iterLHS->first < iterRHS->first)
			{
				cRet.insert(*iterLHS);
				iterLHS++;
			}
			else if (iterRHS->first < iterLHS->first)
			{
				cRet.insert(*iterRHS);
				iterRHS++;
			}
			else
			{
				int iCmpStatus = CompareSchedCond(*iterLHS, *iterRHS);
				if (iCmpStatus == 0)
				{
					cRet.insert(*iterLHS);
				}
				else if (iCmpStatus == -1)
				{
					//LHS stricter	
					cRet.insert(*iterLHS);
					iRelaxed = 0;
				}
				else if (iCmpStatus == 1)
				{
					//LHS relaxed	
					cRet.insert(*iterRHS);
					iStricter = 0;
				}
				iterLHS++;
				iterRHS++;
			}
		}

		if (iterLHS == rcLHS.end())
		{
			//RHS has more conditions.
			for (; iterRHS != rcRHS.end(); iterRHS++)
			{
				if (!IsSchedCondRedundant(*iterRHS))
				{
					cRet.insert(*iterRHS);										
				}
			}
		}

		if (iterRHS == rcRHS.end())
		{
			//LHS has more conditions.
			for (; iterLHS != rcLHS.end(); iterLHS++)
			{
				if (!IsSchedCondRedundant(*iterLHS))
				{					
					cRet.insert(*iterLHS);
				}
			}
		}
		return cRet;
	}

	SchedConds DisjunctSchedConds(SchedConds & rcLHS, SchedConds & rcRHS)
	{
		SchedConds cRet;
		typename SchedConds::iterator iterLHS = rcLHS.begin();
		typename SchedConds::iterator iterRHS = rcRHS.begin();
		int iStricter = 1;
		int iRelaxed = 1;
		for (; iterLHS != rcLHS.end() && iterRHS != rcRHS.end();)
		{
			if (iterLHS->first < iterRHS->first)
			{				
				iterLHS++;
			}
			else if (iterRHS->first < iterLHS->first)
			{				
				iterRHS++;
			}
			else
			{
				int iCmpStatus = CompareSchedCond(*iterLHS, *iterRHS);
				if (iCmpStatus == 0)
				{
					cRet.insert(*iterLHS);
				}
				else if (iCmpStatus == -1)
				{
					//LHS stricter	
					cRet.insert(*iterRHS);
					iRelaxed = 0;
				}
				else if (iCmpStatus == 1)
				{
					//LHS relaxed	
					cRet.insert(*iterLHS);
					iStricter = 0;
				}
				iterLHS++;
				iterRHS++;
			}
		}	
		return cRet;
	}
	
	void setEnforcedOrder(int iTaskA, int iTaskB, int iPartialOrder)
	{
		if (iPartialOrder == -1)
		{
			m_cEnforcedPartialOrder[iTaskA][iTaskB] = -1;
			m_cEnforcedPartialOrder[iTaskB][iTaskA] = -1;
		}
		else
		{
			m_cEnforcedPartialOrder[iTaskA][iTaskB] = iPartialOrder;
			m_cEnforcedPartialOrder[iTaskB][iTaskA] = 1 - iPartialOrder;
		}		
	}

protected:

	void InitiailizePartialOrderTable()
	{
		int iTaskNum = m_pcTaskSet->getNumTasks();
		m_cEnforcedPartialOrder.reserve(iTaskNum);
		for (int i = 0; i < iTaskNum; i++)
		{
			m_cEnforcedPartialOrder.push_back(IntArray1D());
			IntArray1D & rcArray = m_cEnforcedPartialOrder.back();
			rcArray.reserve(iTaskNum);
			for (int j = 0; j < iTaskNum; j++)
			{
				rcArray.push_back(-1);
			}
		}
	}



	virtual bool SchedTest(int iTaskIndex, PriorityAssignment & rcPriorityAssignment, SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed)
	{
		cout << "Base class doesn't assume any system. Don't use!" << endl;
		return false;
	}

	virtual int PreSchedTest(int iTaskIndex, PriorityAssignment & rcPriorityAssignment, bool & rbSchedStatus, SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed)
	{
		return 1;
	}

	virtual int PostSchedTest(int iTaskIndex, PriorityAssignment & rcPriorityAssignment, bool bSchedStatus, SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed)
	{
		return 1;
	}

	virtual void PreComputeUnschedCores(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedConstFixed, SchedCondsArray & rcUnschedCores, int iLimit, double dTimeout = 1e74)
	{
	}

	virtual int ComputeAllCores(SchedConds & rcSchedConst, UnschedCores & rcUnschedCores, int iLimit, double dTimeout = 1e74)
	{

		return 0;
	}

	inline virtual bool IsDistinguishable(SchedCondContent & rcLHS, SchedCondContent & rcRHS)
	{
		return abs(rcRHS - rcLHS) > 1;
	}

	virtual void ConvertToMUS_Schedulable(SchedCond & rcErased, SchedConds & rcSchedCondsFlex, SchedConds & rcShcedCondsFixed, PriorityAssignment & rcPriorityAssignment)
	{
		//rcSchedCondsFlex.insert(rcErased);
		//use bisection method.
		Monotonicity cMono = SchedCondMonotonicity(rcErased.first);
		SchedCondContent cLB = SchedCondLowerBound(rcErased.first);
		SchedCondContent cUB = SchedCondUpperBound(rcErased.first);		
		SchedCondContent & rcUnschedBound = rcErased.second;
		SchedCondContent &  rcSchedBound = (cMono == Up) ? cUB : cLB;
		PriorityAssignment cPriorityAssignmentThis(*m_pcTaskSet);
		PriorityAssignment & cPriorityAssignmentLastUnsched = rcPriorityAssignment;
		cPriorityAssignmentThis.CopyFrom_Strict(cPriorityAssignmentLastUnsched);
		while (IsDistinguishable(rcSchedBound, rcUnschedBound))
		{
			SchedCondContent cMid = (rcUnschedBound + rcSchedBound) / 2;
			rcSchedCondsFlex[rcErased.first] = cMid;
			bool bStatus = IsSchedulable(rcSchedCondsFlex, rcShcedCondsFixed, cPriorityAssignmentThis);
			if (bStatus)
			{
				rcSchedBound = cMid;
				cPriorityAssignmentThis.CopyFrom_Strict(cPriorityAssignmentLastUnsched);
			}
			else
			{
				rcUnschedBound = cMid;
				cPriorityAssignmentLastUnsched.CopyFrom_Strict(cPriorityAssignmentThis);
			}
		}
		rcSchedCondsFlex[rcErased.first] = rcUnschedBound;
	}

	virtual void MUSConvertOrder(SchedConds & rcSchedCondsFlex, vector<SchedCond> & rvectorSchedConds)
	{
		for (typename SchedConds::iterator iter = rcSchedCondsFlex.begin(); iter != rcSchedCondsFlex.end(); iter++)
		{
			rvectorSchedConds.push_back(*iter);
		}
	}

	virtual void PreConvertToMUS(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed, SchedConds & rcCore)
	{

	}

	virtual void PostConvertToMUS(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed, SchedConds & rcCore)
	{

	}

	virtual int ConvertToMUS(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed, SchedConds & rcCore)
	{
		PreConvertToMUS(rcSchedCondsFlex, rcSchedCondsFixed, rcCore);
		if (IsSchedulable(rcSchedCondsFlex, rcSchedCondsFixed))
		{
			return 0;
		}
		rcCore = rcSchedCondsFlex;
		vector<SchedCond> vectorSorted;
		vectorSorted.reserve(rcSchedCondsFlex.size());
		MUSConvertOrder(rcSchedCondsFlex, vectorSorted);
		PriorityAssignment cPriorityAssignment(*m_pcTaskSet);
		PriorityAssignment cPriorityAssignmentLastUnsched(*m_pcTaskSet);
		for (typename vector<SchedCond>::iterator iter = vectorSorted.begin(); iter != vectorSorted.end(); iter++)
		{
			SchedCond cSchedCond = *iter;
			rcCore.erase(iter->first);
			bool bSchedulable = IsSchedulable(rcCore, rcSchedCondsFixed, cPriorityAssignment);
			if (bSchedulable)
			{
				ConvertToMUS_Schedulable(cSchedCond, rcCore, rcSchedCondsFixed, cPriorityAssignmentLastUnsched);
				cPriorityAssignment.CopyFrom_Strict(cPriorityAssignmentLastUnsched);
			}
			else
			{
				cPriorityAssignmentLastUnsched.CopyFrom_Strict(cPriorityAssignment);
			}
		}
		PostConvertToMUS(rcSchedCondsFlex, rcSchedCondsFixed, rcCore);
		return 1;
	}

	virtual bool Recur_MinimalRelaxation(SchedConds & rcSchedCondsFlex, SchedCond cSchedCondInFlex, SchedCond cSchedCondRef)
	{		
		Monotonicity cMono = SchedCondMonotonicity(cSchedCondRef.first);
		if (IsSchedCondRedundant(cSchedCondRef))
			return false;
		SchedCondContent cRelaxed = (cMono == Up) ? max(cSchedCondInFlex.second, cSchedCondRef.second + 1) : min(cSchedCondInFlex.second, cSchedCondRef.second - 1);
		cSchedCondInFlex.second = cRelaxed;
		if (IsSchedCondRedundant(cSchedCondInFlex))
			rcSchedCondsFlex.erase(cSchedCondInFlex.first);
		else
			rcSchedCondsFlex[cSchedCondInFlex.first] = cSchedCondInFlex.second;
		return true;
	}

	int Recur_v2(int iLevel, int iLimit, SchedConds & rcSchedCondsFlex,
		SchedConds & rcSchedCondsFixed, UnschedCores & rdequeCollectedCores, double dTimeStart, double dTimeout)
	{
		if ((double)time(NULL) - dTimeStart >= dTimeout)
			return 0;
		
		if ((iLimit != -1) && (rdequeCollectedCores.size() >= iLimit))
			return 0;

		bool bNewMiniSet = false;
		if (rdequeCollectedCores.size() <= iLevel)
		{
			SchedConds cMUS;
			bNewMiniSet = ConvertToMUS(rcSchedCondsFlex, rcSchedCondsFixed, cMUS);
			if (bNewMiniSet)
			{
				NewCoreCallBack(cMUS);
				rdequeCollectedCores.push_back(cMUS);
#if 0
				cout << "Core " << rdequeCollectedCores.size() << ": ";
				PrintSetContent(rdequeCollectedCores.back().getSetContent());
#endif
			}
			else
			{
				return 0;
			}
		}

		if (bNewMiniSet || rdequeCollectedCores.size() > iLevel)
		{
			SchedConds & rcMUS = rdequeCollectedCores[iLevel];
			//if (!(rcSchedCondsFlex >= rcMUS))
			if (CompareSchedConds(rcSchedCondsFlex, rcMUS) != -1) // SchedConsFlex is not stricter than the core. Meaning it will not leads to rcMUS.
				return -1;
			for (typename SchedConds::iterator iter = rcMUS.begin(); iter != rcMUS.end(); iter++)
			{
				assert(rcSchedCondsFlex.count(iter->first));
				//SchedCondType iSchedCondEle = *rcSchedCondsFlex.find(*iter);
				SchedCond cSchedCondEle = *rcSchedCondsFlex.find(iter->first);
				SchedCond cSchedCondEleMUS = *iter;
				//rcSchedCondsFlex.erase(cSchedCondEle.first);
				if (!Recur_MinimalRelaxation(rcSchedCondsFlex, cSchedCondEle, cSchedCondEleMUS))	continue;
				int iNextLevel = iLevel + 1;
				int iStatus = Recur(iNextLevel, iLimit, rcSchedCondsFlex, rcSchedCondsFixed, rdequeCollectedCores, dTimeStart, dTimeout);
				while (iStatus == -1)
				{
					iNextLevel++;
					iStatus = Recur(iNextLevel, iLimit, rcSchedCondsFlex, rcSchedCondsFixed, rdequeCollectedCores, dTimeStart, dTimeout);
				}
				//rcSchedCondsFlex.insert(cSchedCondEle);
				rcSchedCondsFlex[cSchedCondEle.first] = cSchedCondEle.second;
			}
		}
		return 0;
	}

	int Recur(int iLevel, int iLimit, SchedConds & rcSchedCondsFlex,
		SchedConds & rcSchedCondsFixed, UnschedCores & rdequeCollectedCores, double dTimeStart, double dTimeout)
	{
		if ((double)time(NULL) - dTimeStart >= dTimeout)
			return 0;
		
		if ((iLimit != -1) && (rdequeCollectedCores.size() >= iLimit))
			return 0;

		bool bNewMiniSet = false;
		if ((rdequeCollectedCores.size() == iLevel) || (rdequeCollectedCores.empty()))
		{
			SchedConds cMUS;
			bNewMiniSet = ConvertToMUS(rcSchedCondsFlex, rcSchedCondsFixed, cMUS);
			if (bNewMiniSet)
			{
				NewCoreCallBack(cMUS);
				rdequeCollectedCores.push_back(cMUS);
#if 0
				cout << "Core " << rdequeCollectedCores.size() << ": ";
				PrintSetContent(rdequeCollectedCores.back().getSetContent());
#endif
			}
			else
			{
				return 0;
			}
		}

		if (bNewMiniSet || rdequeCollectedCores.size() > iLevel)
		{
			SchedConds & rcMUS = rdequeCollectedCores[iLevel];
			//if (!(rcSchedCondsFlex >= rcMUS))
			int iNextLevel = iLevel;
			int iCmpStatus = CompareSchedConds(rcSchedCondsFlex, rcMUS);			
			while ((iCmpStatus != -1 ) && (iCmpStatus != 0))//As long as the mus is not stricter or equal 
			{	
				iNextLevel += 1;
				if (iNextLevel == rdequeCollectedCores.size()) break;
				rcMUS = rdequeCollectedCores[iNextLevel];
				iCmpStatus = CompareSchedConds(rcSchedCondsFlex, rcMUS);
			}

			if (iNextLevel < rdequeCollectedCores.size())
			{
				for (typename SchedConds::iterator iter = rcMUS.begin(); iter != rcMUS.end(); iter++)
				{
					assert(rcSchedCondsFlex.count(iter->first));
					//SchedCondType iSchedCondEle = *rcSchedCondsFlex.find(*iter);
					SchedCond cSchedCondEle = *rcSchedCondsFlex.find(iter->first);
					SchedCond cSchedCondEleMUS = *iter;
					//rcSchedCondsFlex.erase(cSchedCondEle.first);
					if (!Recur_MinimalRelaxation(rcSchedCondsFlex, cSchedCondEle, cSchedCondEleMUS))	continue;
					iNextLevel += 1;
					int iStatus = Recur(iNextLevel, iLimit, rcSchedCondsFlex, rcSchedCondsFixed, rdequeCollectedCores, dTimeStart, dTimeout);
					//rcSchedCondsFlex.insert(cSchedCondEle);
					rcSchedCondsFlex[cSchedCondEle.first] = cSchedCondEle.second;
				}
			}
			else
			{
				int iStatus = Recur(iNextLevel, iLimit, rcSchedCondsFlex, rcSchedCondsFixed, rdequeCollectedCores, dTimeStart, dTimeout);
			}

			
		}
		return 0;
	}

	virtual void IsSchedulablePreWork(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed, PriorityAssignment & rcPriorityAssignment)
	{

	}

	virtual int NewCoreCallBack(SchedConds & cMUS)
	{
		return 0;
	}

	virtual void HighestPriorityAssignment_Prework(int iTaskIndex, SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed, PriorityAssignment & rcPriorityAssignment)
	{

	}	

	friend ostream & operator << (ostream & rcOS, SchedCond & rcLHS)
	{
		rcOS << "(";
		rcOS << rcLHS.first;
		rcOS << " : ";
		rcOS << rcLHS.second;
		rcOS << ")";
	}

	friend istream & operator >> (istream & rcIS, SchedCond & rcRHS)
	{
		char xChar = 0;
		rcIS >> xChar;
		if (xChar != '(')
		{
			cout << "Instream General SchedCond error. Do not proceed." << endl;
			while (1);
			return rcIS;
		}
		rcIS >> rcRHS.first;
		rcIS >> xChar;
		if (xChar != ':')
		{
			cout << "Instream General SchedCond error. Do not proceed." << endl;
			while (1);
			return rcIS;
		}
		rcIS >> rcRHS.second;
		rcIS >> xChar;
		if (xChar != ')')	return rcIS;
		return rcIS;
	}
};