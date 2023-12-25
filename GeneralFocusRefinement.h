#pragma once
#include <iostream>
using namespace std;
#include <ilcplex/ilocplex.h>
#include "GeneralUnschedCoreComputer.hpp"
#include "Util.h"
#include "StatisticSet.h"
template < class UnschedCoreComputer_MonoInt_T >
class GeneralFocusRefinement_MonoInt
{
public:
	typedef UnschedCoreComputer_MonoInt_T UnschedCoreComputer_MonoInt;
	typedef typename UnschedCoreComputer_MonoInt::UnschedCore UnschedCore;
	typedef typename UnschedCoreComputer_MonoInt::SchedConds SchedConds;
	typedef typename UnschedCoreComputer_MonoInt::SchedCond SchedCond;
	typedef typename UnschedCoreComputer_MonoInt::SchedCondType SchedCondType;
	typedef typename UnschedCoreComputer_MonoInt::SchedCondContent SchedCondContent;
	typedef typename UnschedCoreComputer_MonoInt::SystemType SystemType;	
	typedef deque<UnschedCore> UnschedCores;


	class SchedCondComparator
	{
	public:
		bool operator() (const SchedCond & rcLHS, const SchedCond & rcRHS)
		{
			if (rcLHS.first < rcRHS.first)
			{
				return true;
			}
			else if (rcLHS.first > rcRHS.first)
			{
				return false;
			}
			else
			{
				if (rcLHS.second < rcRHS.second)
				{
					return true;
				}
				else if (rcLHS.second > rcRHS.second)
				{
					return false;
				}
				else
				{
					return false;
				}
			}
		}
	};

	typedef map < SchedCondType, map<SchedCondContent, IloNumVar> > SchedCondVars;
	typedef map<SchedCondType, IloNumVar> SchedCondTypeVars;

	enum IterationState
	{
		Search,
		Optimal,
		Abort,
		Infeasible,
		Timeout,
		Feasible,
		
	};

	struct IterationStatus
	{
		int iIterationN;
		SchedConds cSchedCondsConfig;
		UnschedCores cThisUnschedCores;
		int iTotalCores;
		double dObjective;
		double dBestFeasible;
		double dObjLB;
		double dTime;
		double dCPUTime;
		double dUCComputeTime;
		double dILPTime;
		IterationState enumState;
		IterationStatus()
		{
			iIterationN = 0;
			iTotalCores = 0;
			dObjective = 0;
			dBestFeasible = 0;
			dTime = 0;
			dCPUTime = 0;
			dUCComputeTime = 0;
			dILPTime = 0;
			enumState = IterationState::Search;
		}
	};

	enum enumAlgConfig
	{
		Lazy,
		Eager
	};

protected:
	SystemType * m_pcTaskSet;
	UnschedCoreComputer_MonoInt m_cUnschedCoreComputer;

	SchedConds m_cSchedCondsConfig;
	deque<IterationStatus> m_dequeIterationLogs;
	string m_stringLogFile;
	double m_dTotalTime;
	PriorityAssignment m_cPriorityAssignmentSol;
	IterationStatus m_cFinalStatistic;
	double m_dStartSeconds;
	double m_dObjective;
	double m_dSolveCPUTime;
	double m_dSolveWallTime;

	//Configuration
	enumAlgConfig m_enumAlgConfig;
	int m_iAlgConfig;
	int m_iSubILPDisplay;
	int m_iCorePerIter;
	double m_dILPGap;	


	//Run Time Data
	double m_dObjUB;
	double m_dObjLB;
	SchedConds m_cUBSchedCondsSolution;
	int m_iBestSolFeasible;
	UnschedCores m_cUnschedCores;
	SchedConds m_cSchedCondsFixed;
	typedef deque< std::pair<double, SchedConds> > SolutionSchedCondsDeque;
	SolutionSchedCondsDeque m_dequeSolutionSchedConds;
	IloCplex::Status m_cCplexStatus;

	//Reference SchedConds Guided Search
	SchedConds m_cReferenceSchedConds;

public:
	GeneralFocusRefinement_MonoInt()
	{
		GenerateFixedSchedConds(m_cSchedCondsFixed);
		m_enumAlgConfig = enumAlgConfig::Lazy;
		m_iSubILPDisplay = 0;
		m_iCorePerIter = 10;
		m_dObjUB = 1e74;
		m_dObjLB = -1e74;
	}

	GeneralFocusRefinement_MonoInt(SystemType & rcTaskSet)
		:m_cUnschedCoreComputer(rcTaskSet)
	{
		m_pcTaskSet = &rcTaskSet;
		GenerateFixedSchedConds(m_cSchedCondsFixed);
		m_enumAlgConfig = enumAlgConfig::Lazy;
		m_iSubILPDisplay = 0;
		m_iCorePerIter = 10;
		m_dILPGap = 1e-4;
		m_dObjUB = 1e74;
		m_dObjLB = -1e74;
	}

	~GeneralFocusRefinement_MonoInt()
	{

	}
	
	int FocusRefinement(int iDisplay, double dTimeout = 1e74)
	{
		assert(m_pcTaskSet);		
		m_cSchedCondsConfig.clear();
		GenerateFixedSchedConds(m_cSchedCondsFixed);
		char axBuffer[256] = { 0 };
		m_dObjUB = getObjUB();
		//m_dObjLB = -1e74;
		m_dStartSeconds = getWallTimeSecond();
		int iStatus = 0;
		//Test if schedulable first	
		if (iDisplay == 2)
			cout << "Testing System Schedulability...." << endl;
		if (!IsSchedulable(m_cSchedCondsConfig, m_cSchedCondsFixed))
		{
			IterationStatus tagStatus;
			tagStatus.dTime = getWallTimeSecond() - m_dStartSeconds;
			tagStatus.enumState = IterationState::Infeasible;
			m_dequeIterationLogs.push_back(tagStatus);
			if (iDisplay)	cout << "Infeasible" << endl;
			PrintIterStatusToLog(m_dequeIterationLogs.size() - 1);
			ComputeStatistic();
			return 0;
		}
		if (iDisplay == 1)
			cout << "system schedulable. start optimization..." << endl;
		Timer cTimer;
		cTimer.Start();
		double dLastObjective = -1e74;
		int iIterationCount = 0;
		double dUCComputeTime = 0;
		bool bCanExit = false;
		int iLastCoreSize = 0;
		while (1)
		{
			bCanExit = false;
			if (iDisplay == 2)	cout << "ILP query...." << endl;
			iStatus = SolveFR(m_iSubILPDisplay, m_cSchedCondsFixed, m_cUnschedCores, m_dObjLB, 1e74, 1e74);
			dLastObjective = m_dObjective;
			double dUCSeconds = getWallTimeSecond();
			UnschedCores cThisCores;
			if (iDisplay == 2)
			{
				cout << "SchedConds: ";
				m_cUnschedCoreComputer.PrintSchedConds(m_cSchedCondsConfig, cout);
				cout << endl;
			}
			if (iDisplay == 2) cout << "Unched-core computation..." << endl;
			switch (m_enumAlgConfig)
			{
			case enumAlgConfig::Lazy:
				ComputeUnschedCores(m_cSchedCondsConfig, m_cSchedCondsFixed, cThisCores, m_iCorePerIter);
				break;
			case enumAlgConfig::Eager:
				ComputeUnschedCoreForAllSol(m_cSchedCondsFixed, cThisCores, m_iCorePerIter);
				break;
			default:
				cout << "Unknown Algorithm Configuration" << endl;
				assert(0);
				break;
			}
			if (iDisplay == 2) cout << "Unched-core computation complete...." << endl;
			PostComputeUnschedCore(cThisCores);
			for (typename UnschedCores::iterator iter = cThisCores.begin();
				iter != cThisCores.end(); iter++)
			{
				m_cUnschedCores.push_back(*iter);
			}

			//Update Iteration Status Statistic
			double dUCComputeTime = (double)(getWallTimeSecond() - dUCSeconds);
			double dTimeByNow = (double)(getWallTimeSecond() - m_dStartSeconds);
			IterationStatus tagStatus;
			iIterationCount++;
			tagStatus.cSchedCondsConfig = m_cSchedCondsConfig;
			tagStatus.iIterationN = iIterationCount;
			tagStatus.dObjective = m_dObjective;
			tagStatus.dBestFeasible = m_dObjUB;
			tagStatus.dObjLB = m_dObjLB;
			tagStatus.iTotalCores = m_cUnschedCores.size();
			tagStatus.cThisUnschedCores = cThisCores;
			tagStatus.dTime = dTimeByNow;
			tagStatus.dILPTime = m_dSolveWallTime;
			tagStatus.dUCComputeTime = dUCComputeTime;
			cTimer.Stop();
			tagStatus.dCPUTime = cTimer.getElapsedTime_ms();
			tagStatus.enumState = IterationState::Search;
			m_dequeIterationLogs.push_back(tagStatus);
			//Handle the result					
			if (m_iBestSolFeasible == 1)
			{
				//A feasible solution found
				if (iDisplay == 2)	cout << "optimal solution found" << endl;
				m_cPriorityAssignmentSol = PriorityAssignment(*m_pcTaskSet);
				IsSchedulable(m_cSchedCondsConfig, m_cSchedCondsFixed, m_cPriorityAssignmentSol);
				if ((m_cCplexStatus == IloCplex().Optimal) || (m_cCplexStatus == IloCplex().OptimalTol))
				{
					m_dequeIterationLogs.back().enumState = Optimal;
					bCanExit = true;
					iStatus = 1;
				}
				else
				{
					cout << m_cCplexStatus << endl;
					cout << "Best feasible yet Cplex not indicating optimality?" << endl;
					assert(0);
				}

				BestFeasibleHook(bCanExit);
			}
			else if (iStatus == 0)
			{
				m_dequeIterationLogs.back().enumState = IterationState::Infeasible;
				bCanExit = true;
				iStatus = 0;
			}
			else if (dTimeByNow > dTimeout)
			{
				m_dequeIterationLogs.back().enumState = Timeout;
				bCanExit = true;
				iStatus = -1;
			}

			//Display Status
			if (iDisplay)
			{
				PrintIterationStatusSimple(m_dequeIterationLogs.back());
				//system("pause");
			}

			PrintIterStatusToLog(m_dequeIterationLogs.size() - 1);
			if (bCanExit)
			{
				if (iDisplay) cout << endl;
				break;
			}


		}

		ComputeStatistic();
		if (iStatus == 1)
		{
			OptimalHook();
		}
		return iStatus;
	}

	void setLogFile(char axCommand[])
	{
		m_stringLogFile = string(axCommand);
		if (!m_stringLogFile.empty())
		{
			ofstream ofstreamLogFile(m_stringLogFile.data(), ios::out);
			ofstreamLogFile.close();
		}
	}

	PriorityAssignment getSolPA()
	{
		return m_cPriorityAssignmentSol;
	}

	StatisticSet GenerateStatisticFile(const char axFileName[])
	{
		StatisticSet cStatistic;
		cStatistic.setItem("Total Time (ms)", m_cFinalStatistic.dTime);
		cStatistic.setItem("Total CPU Time (ms)", m_cFinalStatistic.dCPUTime);
		cStatistic.setItem("Iteration", m_cFinalStatistic.iIterationN);
		cStatistic.setItem("Cores", m_cFinalStatistic.iTotalCores);
		cStatistic.setItem("Total UC Time (ms)", m_cFinalStatistic.dUCComputeTime);
		cStatistic.setItem("Total ILP Time (ms)", m_cFinalStatistic.dILPTime);
		cStatistic.setItem("Objective", m_dObjective);
		cStatistic.setItem("Best Feasible", m_cFinalStatistic.dBestFeasible);
		//cStatistic.setItem("Objective UB", m_pcTaskSet->getWorstDelay());
		cStatistic.setItem("Objective UB", getObjUB());

		switch (m_cFinalStatistic.enumState)
		{
		case IterationState::Search:
			cStatistic.setItem("Status", "Searching");
			break;
		case IterationState::Abort:
			cStatistic.setItem("Status", "Abort");
			break;
		case IterationState::Optimal:
			cStatistic.setItem("Status", "Optimal");
			break;
		case IterationState::Infeasible:
			cStatistic.setItem("Status", "Infeasible");
			break;
		case IterationState::Timeout:
			cStatistic.setItem("Status", "Timeout");
			break;
		default:
			my_assert(0, "Unknown Status");
			break;
		}

		//cStatistic.setItem("", m_cFinalStatistic);
		if (strlen(axFileName) == 0)
			return cStatistic;
		char axBuffer[256] = { 0 };
		sprintf(axBuffer, "%s.txt", axFileName);
		cStatistic.WriteStatisticText(axBuffer);

		sprintf(axBuffer, "%s.rslt", axFileName);
		cStatistic.WriteStatisticImage(axBuffer);
		return cStatistic;
	}

	bool IsSchedulable(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed)
	{
		return m_cUnschedCoreComputer.IsSchedulable(rcSchedCondsFlex, rcSchedCondsFixed);
	}

	bool IsSchedulable(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed, PriorityAssignment & rcPriorityAssignment)
	{
		return m_cUnschedCoreComputer.IsSchedulable(rcSchedCondsFlex, rcSchedCondsFixed, rcPriorityAssignment);
	}

	virtual int ComputeUnschedCores(SchedConds & rcSchedCondsFlex, SchedConds & rcSchedCondsFixed, UnschedCores & rcUnschedCores, int iLimit)
	{
		int iNum = m_cUnschedCoreComputer.ComputeUnschedCores(rcSchedCondsFlex, rcSchedCondsFixed, rcUnschedCores, iLimit, 1e74);
		return iNum;
	}

	void setAlgorithm(int iAlgConfig)	{ m_iAlgConfig = iAlgConfig; }

	void setAlgorithm(enumAlgConfig iAlgConfig)	{ m_enumAlgConfig = iAlgConfig; }

	void setSubILPDisplay(int iDisplay)	{ m_iSubILPDisplay = iDisplay; }

	void setCorePerIter(int iCorePerIter)	{ m_iCorePerIter = iCorePerIter; }

	void setILPGap(double dGap)	{ m_dILPGap = dGap; }

	void setLB(double dLB)	{ m_dObjLB = dLB; }

	void WriteUnschedCoreToFile(char axFileName[])
	{

	}

	void ReadUnschedCoreFromFile(const char axFileName[], UnschedCores & rcUCs)
	{
		ifstream ifstreamInFile(axFileName, ios::in);
		char axBuffer[8192] = { 0 };
		char axCoreKey[] = "Core:";
		while (ifstreamInFile.eof() == false)
		{
			ifstreamInFile.getline(axBuffer, 8192);
			memcpy(axCoreKey, axBuffer, 5);
			if (strcmp(axCoreKey, "Core:") != 0)	continue;
			SchedConds cNewUC;
			m_cUnschedCoreComputer.ReadSchedCondsFromStr(cNewUC, &axBuffer[5]);
			rcUCs.push_back(cNewUC);
		}		
	}

	PriorityAssignment & getPriorityAssignmentSol()	{ return m_cPriorityAssignmentSol; }

	virtual UnschedCoreComputer_MonoInt & getUnschedCoreComputer()	{ return m_cUnschedCoreComputer; }

	SchedConds getOptimalSchedCondConfig()	{ return m_cSchedCondsConfig; }

	void AddFixedSchedCond(SchedCond cFixedSchedCond)
	{

	}

	virtual void ClearFR()
	{
		m_cSchedCondsFixed.clear();
		GenerateFixedSchedConds(m_cSchedCondsFixed);
		m_enumAlgConfig = enumAlgConfig::Lazy;
		m_iSubILPDisplay = 0;
		m_iCorePerIter = 10;
		m_dILPGap = 1e-4;
		m_dObjUB = 1e74;
		m_dObjLB = -1e74;

		m_cSchedCondsConfig.clear();
		m_cUnschedCores.clear();
		m_cUBSchedCondsSolution.clear();
		m_cSchedCondsConfig.clear();
		m_cPriorityAssignmentSol = PriorityAssignment(*m_pcTaskSet);
		m_cFinalStatistic = IterationStatus();
		m_dTotalTime = 0;
		m_dequeSolutionSchedConds.clear();
		m_dequeIterationLogs.clear();
		m_iBestSolFeasible = 0;

	}

	double getObj()	{ return m_dObjective; }

	void setReferenceSchedConds(SchedConds & rcSchedConds)	{ m_cReferenceSchedConds = rcSchedConds; }

	void LoadUnschedCores(UnschedCores & rcUnSchedCores)
	{
		m_cUnschedCores.insert(m_cUnschedCores.end(), rcUnSchedCores.begin(), rcUnSchedCores.end());
	}

	const UnschedCores & getUnschedCores()	
	{
		return m_cUnschedCores;
	}

	void getSolutionTrace(deque<SchedConds> & rdequeSolutionTrace)
	{
		for (typename deque<IterationStatus>::iterator iter = m_dequeIterationLogs.begin(); iter != m_dequeIterationLogs.end(); iter++)
		{
			rdequeSolutionTrace.push_back(iter->cSchedCondsConfig);
		}
	}

protected:

	virtual void PostComputeUnschedCore(UnschedCores & rcUnschedCores)
	{
		return;
	}

	void CreateSchedCondVars(const IloEnv & rcEnv, UnschedCores & rcUnschedCores, SchedCondVars & rcSchedCondVars)
	{
		for (typename UnschedCores::iterator iterCore = rcUnschedCores.begin(); iterCore != rcUnschedCores.end(); iterCore++)
		{
			for (typename UnschedCore::iterator iterEle = iterCore->begin(); iterEle != iterCore->end(); iterEle++)
			{
				rcSchedCondVars[iterEle->first][iterEle->second] = IloNumVar(rcEnv, 0.0, 1.0, IloNumVar::Int);
			}
		}
	}

	void CreateSchedCondTypeVars(const IloEnv & rcEnv, UnschedCores & rcUnschedCores, SchedCondTypeVars & rcSchedCondTypeVars)
	{
		for (typename UnschedCores::iterator iterCore = rcUnschedCores.begin(); iterCore != rcUnschedCores.end(); iterCore++)
		{
			for (typename UnschedCore::iterator iterEle = iterCore->begin(); iterEle != iterCore->end(); iterEle++)
			{				
				if (rcSchedCondTypeVars.count(iterEle->first) == 0)
				{
					SchedCondContent cLB = m_cUnschedCoreComputer.SchedCondLowerBound(iterEle->first);
					SchedCondContent cUB = m_cUnschedCoreComputer.SchedCondUpperBound(iterEle->first);
					rcSchedCondTypeVars[iterEle->first] = IloNumVar(rcEnv, cLB, cUB, IloNumVar::Float);
				}
			}
		}
	}

	void GenSchedCondTypeConst(const IloEnv & rcEnv, SchedCondTypeVars & rcSchedCondTypeVars, SchedCondVars & rcSchedCondVars, IloRangeArray & rcConst)
	{
		//typedef map < SchedCondType, map<SchedCondContent, IloNumVar> > SchedCondVars;
		//typedef map<SchedCondType, IloNumVar> SchedCondTypeVars;
		for (typename map < SchedCondType, map<SchedCondContent, IloNumVar> >::iterator iterType = rcSchedCondVars.begin();
			iterType != rcSchedCondVars.end(); iterType++)
		{
			typename UnschedCoreComputer_MonoInt::Monotonicity cMono = m_cUnschedCoreComputer.SchedCondMonotonicity(iterType->first);
			assert(rcSchedCondTypeVars.count(iterType->first));
			for (typename map < SchedCondContent, IloNumVar >::iterator iterContent = iterType->second.begin();
				iterContent != iterType->second.end(); iterContent++)
			{				
				if (cMono == UnschedCoreComputer_MonoInt::Up)
				{					
					rcConst.add(IloRange(rcEnv,
						0.0,
						rcSchedCondTypeVars[iterType->first] - (iterContent->first + 1) * (1.0 - iterContent->second),
						IloInfinity));
				}
				else
				{
					SchedCondContent cUB = m_cUnschedCoreComputer.SchedCondUpperBound(iterType->first);
					rcConst.add(IloRange(rcEnv, 
						-IloInfinity, 
						rcSchedCondTypeVars[iterType->first] - (iterContent->first - 1) - cUB * (iterContent->second),
						0.0));
				}
			}
		}
	}

	void GenSchedCondValueRelationalConst(const IloEnv & rcEnv, SchedCondVars & rcSchedCondVars, IloRangeArray & rcConst)
	{
		for (typename map < SchedCondType, map<SchedCondContent, IloNumVar> >::iterator iterType = rcSchedCondVars.begin();
			iterType != rcSchedCondVars.end(); iterType++)
		{
			typename UnschedCoreComputer_MonoInt::Monotonicity cMono = m_cUnschedCoreComputer.SchedCondMonotonicity(iterType->first);			
			typename map < SchedCondContent, IloNumVar >::iterator iterContent = iterType->second.begin();
			typename map < SchedCondContent, IloNumVar >::iterator iterContentNext = iterContent;
			iterContentNext++;
			for (; iterContentNext != iterType->second.end(); iterContent++, iterContentNext++)
			{
				if (cMono == UnschedCoreComputer_MonoInt::Up)
				{
					rcConst.add(IloRange(rcEnv,
						-IloInfinity,
						iterContent->second - iterContentNext->second,
						0.0));
				}
				else
				{					
					rcConst.add(IloRange(rcEnv,
						0.0,
						iterContent->second - iterContentNext->second,
						IloInfinity));
				}
			}
		}
	}

	void GenUnschedCoreConst(const IloEnv & rcEnv, UnschedCores & rcUnschedCores, SchedCondVars & rcSchedCondVars, IloRangeArray & rcConst)
	{
		for (typename UnschedCores::iterator iterCore = rcUnschedCores.begin(); iterCore != rcUnschedCores.end(); iterCore++)
		{
			IloExpr cExpr(rcEnv);
			cExpr += 0;
			for (typename UnschedCore::iterator iterEle = iterCore->begin(); iterEle != iterCore->end(); iterEle++)
			{
				assert(rcSchedCondVars.count(iterEle->first));
				assert(rcSchedCondVars[iterEle->first].count(iterEle->second));
				cExpr += rcSchedCondVars[iterEle->first][iterEle->second];
			}
			rcConst.add(cExpr <= (int)iterCore->size() - 1);
			//rcConst.add(IloRange(rcEnv, 0.0, cExpr, (int)iterCore->size() - 1));						
		}
	}

	IloExpr ScheCondsSatSumExpr(const IloEnv & rcEnv, SchedConds & rcSchedConds, SchedCondVars & rcSchedCondVars)
	{
		IloExpr cSatSum(rcEnv);
		for (typename SchedConds::iterator iter = rcSchedConds.begin(); iter != rcSchedConds.end(); iter++)
		{
			assert(rcSchedCondVars.count(*iter));
			cSatSum += rcSchedCondVars[*iter];
		}
		return cSatSum;
	}

	virtual void GenerateFixedSchedConds(SchedConds & rcSchedCondFixed)
	{
		rcSchedCondFixed = m_cSchedCondsFixed;
		return;
	}

	virtual double getObjUB()
	{
		return 1e74;
	}

	virtual void OptimalHook()
	{

	}

	virtual void setSolverParam(IloCplex & rcSolver)
	{
		rcSolver.setParam(IloCplex::Param::Emphasis::MIP, 0);
		rcSolver.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 1e-9);
		rcSolver.setParam(IloCplex::Param::MIP::Tolerances::Integrality, 0);
	}

	void EnableSolverDisp(int iLevel, IloCplex & rcSolver)
	{
		const IloEnv & rcEnv = rcSolver.getEnv();
		if (iLevel)
		{
			rcSolver.setParam(IloCplex::Param::MIP::Display, iLevel);
		}
		else
		{
			rcSolver.setParam(IloCplex::Param::MIP::Display, 0);
			rcSolver.setWarning(rcEnv.getNullStream());
			rcSolver.setOut(rcEnv.getNullStream());
		}
	}

	virtual int SolveFR(int iDisplay, SchedConds & rcSchedCondFixed, UnschedCores & rcGenUC, double dObjLB, double dObjUB, double dTimeout = 1e74)
	{
		try
		{
			IloEnv cSolveEnv;
			IloModel cModel(cSolveEnv);
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
			if (bStatus)
			{

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

	void ExtractSolutiontPrimitive(IloCplex & rcSolver, int iSolIndex, SchedCondVars & rcSchedCondVars, SchedConds & rcSchedCondsConfig)
	{
		for (typename SchedCondVars::iterator iter = rcSchedCondVars.begin(); iter != rcSchedCondVars.end(); iter++)
		{
			int iValue = round(rcSolver.getValue(iter->second));
			if (iValue == 1)
			{
				rcSchedCondsConfig.insert(iter->first);
			}
		}
	}

	virtual void ExtractSolution(IloCplex & rcSolver, int iSolIndex, SchedCondTypeVars & rcSchedCondTypeVars, SchedConds & rcSchedCondsConfig)
	{		
		for (typename SchedCondTypeVars::iterator iter = rcSchedCondTypeVars.begin(); iter != rcSchedCondTypeVars.end(); iter++)
		{
			double dValue = rcSolver.getValue(iter->second);
			dValue = round(dValue);
			if (m_cUnschedCoreComputer.IsSchedCondRedundant({ iter->first, dValue })) continue;
			rcSchedCondsConfig[iter->first] = dValue;
		}
	}

	virtual void ExtractAllSolution(IloCplex & rcSolver, SchedCondTypeVars & rcSchedCondTypeVars)
	{
		int iSolNum = rcSolver.getSolnPoolNsolns();
		m_dequeSolutionSchedConds.clear();
		for (int i = 0; i < iSolNum; i++)
		{
			SchedConds rcConfig;
			ExtractSolution(rcSolver, i, rcSchedCondTypeVars, rcConfig);
			double dObj = rcSolver.getObjValue(i);
			m_dequeSolutionSchedConds.push_back({ dObj, rcConfig });
		}
	}

	void MarshalAllSolution(SolutionSchedCondsDeque & rdequeSolutionSchedConds)
	{
		if (rdequeSolutionSchedConds.empty())
		{
			cout << "No solution? Is it a linear programming problem?" << endl;
			assert(0);
		}
		else
		{
			assert(!m_dequeSolutionSchedConds.empty());
			m_cSchedCondsConfig = m_dequeSolutionSchedConds.front().second;
		}

		for (typename SolutionSchedCondsDeque::iterator iter = rdequeSolutionSchedConds.begin();
			iter != rdequeSolutionSchedConds.end(); iter++)
		{
			SchedConds & cSchedCondsConfig = iter->second;
			double m_dThisObj = iter->first;
			if (m_dThisObj < m_dObjUB)
			{
				if (IsSchedulable(cSchedCondsConfig, m_cSchedCondsFixed))
				{
					if (iter == rdequeSolutionSchedConds.begin())
					{
						m_iBestSolFeasible = 1;
					}
					m_dObjUB = m_dThisObj;
					m_cUBSchedCondsSolution = cSchedCondsConfig;
				}
			}
			else
				break;
		}

		m_dObjective = rdequeSolutionSchedConds.begin()->first;
		if (m_dObjective > m_dObjLB)
			m_dObjLB = m_dObjective;
		if ((abs(m_dObjUB - m_dObjective) < 1e-5))
		{
			m_cSchedCondsConfig = m_cUBSchedCondsSolution;
			m_iBestSolFeasible = 1;
			return;
		}

		if ((m_iBestSolFeasible == 0 && IsSchedulable(m_cSchedCondsConfig, m_cSchedCondsFixed)))
		{
			m_iBestSolFeasible = 1;
		}
	}

	void ComputeUnschedCoreForAllSol(SchedConds & rcSchedCondsFixed, UnschedCores & rcUnscheCores, int iLimit)
	{
		//int iSize = min(m_dequeSolutionSchedConds.size(), 100);
		int iSize = m_dequeSolutionSchedConds.size();
		for (int i = 0; i < iSize; i++)
		{
			SchedConds & rcSchedCondsConfig = m_dequeSolutionSchedConds[i].second;
			UnschedCores cThisCores;
			ComputeUnschedCores(rcSchedCondsConfig, rcSchedCondsFixed, cThisCores, iLimit);
			int iCoreSize = cThisCores.size();
			for (int j = 0; j < iCoreSize; j++)
			{
				rcUnscheCores.push_back(cThisCores[j]);
			}
		}
	}

	void PrintIterationStatus(int iIndex)
	{
		assert(m_dequeIterationLogs.size() > iIndex);
		IterationStatus & rtagStatus = m_dequeIterationLogs[iIndex];
		cout << "----------------------------------------------" << endl;
		cout << "Iteration " << rtagStatus.iIterationN << endl;
		cout << "Objective: " << rtagStatus.dObjective << endl;
		if (iIndex > 0)
		{
			double dConvergeRate = (rtagStatus.dObjective - m_dequeIterationLogs[iIndex - 1].dObjective) / m_dequeIterationLogs[iIndex - 1].dObjective;
			dConvergeRate = min(dConvergeRate, 1.0);
			printf("Converge Rate: %.2f%%\n", dConvergeRate * 100);
		}

		cout << "Unsched-Cores: " << rtagStatus.iTotalCores << endl;
		int iSize = rtagStatus.cThisUnschedCores.size();
		for (int i = 0; i < iSize; i++)
		{
			cout << "Core: ";
			//PrintSetContent(rtagStatus.cThisUnschedCores[i].getSetContent());
			m_cUnschedCoreComputer.PrintSchedConds(rtagStatus.cThisUnschedCores[i]);
		}
	}

	void PrintIterationStatusSimple(IterationStatus & rcStatus)
	{
#if 0
		cout <<
			"\rIt: " << rcStatus.iIterationN <<
			" | UB: " << rcStatus.dBestFeasible <<
			" | LB: " << rcStatus.dObjLB <<
			" | Obj: " << rcStatus.dObjective <<
			" | UC: " << rcStatus.iTotalCores <<
			" | " << StatusToString(rcStatus.enumState).data() <<
			//" | " << rcStatus.enumState <<
			" | t: " << rcStatus.dTime << "s   ";
#else
		printf("\rIt: %d | UB: %.2e | LB: %.2f | obj: %.2f | UC: %d | %s | t: %.2fs      ",
			rcStatus.iIterationN, rcStatus.dBestFeasible, rcStatus.dObjLB,
			rcStatus.dObjective, rcStatus.iTotalCores, StatusToString(rcStatus.enumState).data(), rcStatus.dTime);
#endif
	}

	void PrintIterationStatus_BinaryFR(int iIndex, int iBinaryIteration, double dObjLB, double dObjUB)
	{
		assert(m_dequeIterationLogs.size() > iIndex);
		IterationStatus & rtagStatus = m_dequeIterationLogs[iIndex];
		cout << "----------------------------------------------" << endl;
		cout << "Binary Iteration: " << iBinaryIteration << endl;
		cout << "Objective Interval: " << "[" << dObjLB << ", " << dObjUB << "]" << endl;;
		cout << "Iteration " << rtagStatus.iIterationN << endl;
		cout << "Objective: " << rtagStatus.dObjective << endl;
		if (iIndex > 0)
		{
			double dConvergeRate = (rtagStatus.dObjective - m_dequeIterationLogs[iIndex - 1].dObjective) / m_dequeIterationLogs[iIndex - 1].dObjective;
			dConvergeRate = min(dConvergeRate, 1.0);
			printf("Converge Rate: %.2f%%\n", dConvergeRate * 100);
		}

		cout << "Unsched-Cores: " << rtagStatus.iTotalCores << endl;
		int iSize = rtagStatus.cThisUnschedCores.size();
		for (int i = 0; i < iSize; i++)
		{
			cout << "Core: ";
			//PrintSetContent(rtagStatus.cThisUnschedCores[i].getSetContent());
			m_cUnschedCoreComputer.PrintSchedConds(rtagStatus.cThisUnschedCores[i]);
		}
	}

	void PrintIteratationStatusSimple_BinaryFR(IterationStatus & rcStatus, int iBinaryIteration, double dObjLB, double dObjUB)
	{
		printf("\rBinStat: %d, [%.2e, %.2e] | It: %d | UB: %.2e | LB: %.2e | obj: %.2e | UC: %d | %s | t: %.2fs      ",
			iBinaryIteration, dObjLB, dObjUB,
			rcStatus.iIterationN, rcStatus.dBestFeasible, rcStatus.dObjLB,
			rcStatus.dObjective, rcStatus.iTotalCores, StatusToString(rcStatus.enumState).data(), rcStatus.dTime);
	}

	void PrintIterStatusToLog_BinaryFR(int iIndex, int iBinaryIteration, double dObjLB, double dObjUB)
	{
		assert(m_dequeIterationLogs.size() > iIndex);
		if (m_stringLogFile.empty())	return;
		char axBuffer[128] = { 0 };
		ofstream ofstreamLogFile(m_stringLogFile.data(), ios::out | ios::app);

		IterationStatus & rtagStatus = m_dequeIterationLogs[iIndex];
		ofstreamLogFile << "----------------------------------------------" << endl;
		ofstreamLogFile << "Binary Iteration: " << iBinaryIteration << endl;
		ofstreamLogFile << "Objective Interval: " << "[" << dObjLB << ", " << dObjUB << "]" << endl;;
		ofstreamLogFile << "Iteration " << rtagStatus.iIterationN << endl;
		ofstreamLogFile << "Objective: " << rtagStatus.dObjective << scientific << endl;
		ofstreamLogFile << "Best Feasible: " << rtagStatus.dBestFeasible << endl;
		ofstreamLogFile << "Objective LB: " << rtagStatus.dObjLB << endl;
		ofstreamLogFile << "Time: " << rtagStatus.dTime << endl;
		ofstreamLogFile << "UC Compute Time: " << rtagStatus.dUCComputeTime << endl;
		ofstreamLogFile << "ILP Time: " << rtagStatus.dILPTime << endl;
		if (iIndex > 0)
		{
			double dConvergeRate = (rtagStatus.dObjective - m_dequeIterationLogs[iIndex - 1].dObjective) / m_dequeIterationLogs[iIndex - 1].dObjective;
			dConvergeRate = min(dConvergeRate, 1.0);
			sprintf(axBuffer, "Converge Rate: %.2f%%\n", dConvergeRate * 100);
			ofstreamLogFile << axBuffer;
		}

		ofstreamLogFile << "SchedCond Configuration: ";
		for (typename UnschedCore::iterator iter = rtagStatus.cSchedCondsConfig.begin();
			iter != rtagStatus.cSchedCondsConfig.end(); iter++)
		{
			m_cUnschedCoreComputer.PrintSchedCond(*iter, ofstreamLogFile);
			ofstreamLogFile << ", ";
			//ofstreamLogFile << *iter << ", ";
		}
		ofstreamLogFile << endl;

		ofstreamLogFile << "Unsched-Cores: " << rtagStatus.iTotalCores << endl;
		int iSize = rtagStatus.cThisUnschedCores.size();
		for (int i = 0; i < iSize; i++)
		{
			ofstreamLogFile << "Core: ";
			for (typename UnschedCore::iterator iter = rtagStatus.cThisUnschedCores[i].begin();
				iter != rtagStatus.cThisUnschedCores[i].end(); iter++)
			{
				//ofstreamLogFile << *iter << ", ";
				m_cUnschedCoreComputer.PrintSchedCond(*iter, ofstreamLogFile);
				ofstreamLogFile << ", ";
			}
			ofstreamLogFile << endl;
		}

		ofstreamLogFile << "Status: ";
		ofstreamLogFile << StatusToString(rtagStatus.enumState);
		ofstreamLogFile << endl;
		PrintExtraInfoToLog(ofstreamLogFile);
		ofstreamLogFile.close();
	}

	void PrintIterStatusToLog(int iIndex)
	{
		assert(m_dequeIterationLogs.size() > iIndex);
		if (m_stringLogFile.empty())	return;
		char axBuffer[128] = { 0 };
		ofstream ofstreamLogFile(m_stringLogFile.data(), ios::out | ios::app);

		IterationStatus & rtagStatus = m_dequeIterationLogs[iIndex];
		ofstreamLogFile << "----------------------------------------------" << endl;
		ofstreamLogFile << "Iteration " << rtagStatus.iIterationN << endl;
		ofstreamLogFile << "Objective: " << rtagStatus.dObjective << scientific << endl;
		ofstreamLogFile << "Best Feasible: " << rtagStatus.dBestFeasible << endl;
		ofstreamLogFile << "Objective LB: " << rtagStatus.dObjLB << endl;
		ofstreamLogFile << "Time: " << rtagStatus.dTime << endl;
		ofstreamLogFile << "UC Compute Time: " << rtagStatus.dUCComputeTime << endl;
		ofstreamLogFile << "ILP Time: " << rtagStatus.dILPTime << endl;
		if (iIndex > 0)
		{
			double dConvergeRate = (rtagStatus.dObjective - m_dequeIterationLogs[iIndex - 1].dObjective) / m_dequeIterationLogs[iIndex - 1].dObjective;
			dConvergeRate = min(dConvergeRate, 1.0);
			sprintf(axBuffer, "Converge Rate: %.2f%%\n", dConvergeRate * 100);
			ofstreamLogFile << axBuffer;
		}

		ofstreamLogFile << "SchedCond Configuration: ";
		for (typename UnschedCore::iterator iter = rtagStatus.cSchedCondsConfig.begin();
			iter != rtagStatus.cSchedCondsConfig.end(); iter++)
		{
			m_cUnschedCoreComputer.PrintSchedCond(*iter, ofstreamLogFile);
			ofstreamLogFile << ", ";
			//ofstreamLogFile << *iter << ", ";
		}
		ofstreamLogFile << endl;

		ofstreamLogFile << "Unsched-Cores: " << rtagStatus.iTotalCores << endl;
		int iSize = rtagStatus.cThisUnschedCores.size();
		for (int i = 0; i < iSize; i++)
		{
			ofstreamLogFile << "Core: ";
			for (typename UnschedCore::iterator iter = rtagStatus.cThisUnschedCores[i].begin();
				iter != rtagStatus.cThisUnschedCores[i].end(); iter++)
			{
				//ofstreamLogFile << *iter << ", ";
				m_cUnschedCoreComputer.PrintSchedCond(*iter, ofstreamLogFile);
				ofstreamLogFile << ", ";
			}
			ofstreamLogFile << endl;
		}

		ofstreamLogFile << "Status: ";
		ofstreamLogFile << StatusToString(rtagStatus.enumState);
		ofstreamLogFile << endl;
		PrintExtraInfoToLog(ofstreamLogFile);
		ofstreamLogFile.close();
	}

	virtual void PrintExtraInfoToLog(ofstream & ofstreamLogFile)
	{

	}

	virtual void ReadLogInfo(const char axFileName[], deque<IterationStatus>  & rdequeIterationStatus)
	{
		ifstream ifstreamInputFile(axFileName, ios::in);
	}

	string StatusToString(IterationState iStatus)
	{
		switch (iStatus)
		{
		case IterationState::Search:
			return string("Searching");
			break;
		case IterationState::Abort:
			return string("Abort");
			break;
		case IterationState::Optimal:
			return string("Optimal");
			break;
		case IterationState::Infeasible:
			return string("Unschedulable");
			break;
		case IterationState::Timeout:
			return string("Timeout");
			break;
		case IterationState::Feasible:
			return string("Feasible");
		default:
			my_assert(0, "Unknown Status");
			break;
		}
		return string("Unknown");
	}

	void ComputeStatistic()
	{
		IterationStatus & rcLast = m_dequeIterationLogs.back();
		m_cFinalStatistic.dObjective = rcLast.dObjective;
		m_cFinalStatistic.iIterationN = rcLast.iIterationN;
		m_cFinalStatistic.cSchedCondsConfig = rcLast.cSchedCondsConfig;
		m_cFinalStatistic.iTotalCores = m_cUnschedCores.size();
		m_cFinalStatistic.dTime = rcLast.dTime * 1000;
		m_cFinalStatistic.dCPUTime = rcLast.dCPUTime * 1000;
		m_cFinalStatistic.enumState = rcLast.enumState;
		m_cFinalStatistic.dBestFeasible = rcLast.dBestFeasible;
		m_cFinalStatistic.dILPTime = 0;
		m_cFinalStatistic.dUCComputeTime = 0;
		for (typename deque<IterationStatus>::iterator iter = m_dequeIterationLogs.begin();
			iter != m_dequeIterationLogs.end(); iter++)
		{
			m_cFinalStatistic.dILPTime += iter->dILPTime;
			m_cFinalStatistic.dUCComputeTime += iter->dUCComputeTime;
		}
		m_cFinalStatistic.dILPTime *= 1000;
		m_cFinalStatistic.dUCComputeTime *= 1000;
	}

	virtual int BinaryFR_BinaryIterationEnd()
	{
		return 0;
	}

	virtual int BinaryFR_RefineBound(double & rdLB, double & rdUB, int iStatus)
	{
		return 0;
	}

	virtual void GenRefSchedCondsModel(const IloEnv & rcEnv, SchedConds & rcRefSchedConds, SchedCondTypeVars & rcSchedCondTypeVars, IloNumVarArray & rcAbsDiffVars, IloRangeArray & rcConst, IloObjective & rcObj)
	{
		rcAbsDiffVars = IloNumVarArray(rcEnv);
		for (typename SchedConds::iterator iter = rcRefSchedConds.begin(); iter != rcRefSchedConds.end(); iter++)
		{
			rcAbsDiffVars.add(IloNumVar(rcEnv, 0.0, IloInfinity));
			IloNumVar & rcThisVar = rcAbsDiffVars[rcAbsDiffVars.getSize() - 1];
			IloNumVar & rcThisSchedCondTypeVar = rcSchedCondTypeVars[iter->first];
			rcConst.add(IloRange(rcEnv,
				0.0, 
				rcThisVar - (rcThisSchedCondTypeVar - iter->second),
				IloInfinity
				));
			rcConst.add(IloRange(rcEnv,
				0.0,
				rcThisVar - (iter->second - rcThisSchedCondTypeVar),
				IloInfinity
				));			
		}
		IloExpr cObjExpr(rcEnv);
		for (int i = 0; i < rcAbsDiffVars.getSize(); i++)
		{
			cObjExpr += rcAbsDiffVars[i];
		}
		rcObj = IloMinimize(rcEnv, cObjExpr);
	}

	virtual void GenRefGuidedSearchModel(const IloEnv & rcEnv, SchedConds & rcRefSchedConds, SchedCondTypeVars & rcSchedCondTypeVars, IloModel & rcModel)
	{
		IloNumVarArray rcAbsDiffVars(rcEnv);
		IloRangeArray rcConst(rcEnv);
		IloObjective rcObj(rcEnv);
		rcAbsDiffVars = IloNumVarArray(rcEnv);
		rcAbsDiffVars.add(IloNumVar(rcEnv, 0.0, IloInfinity));
		for (typename SchedConds::iterator iter = rcRefSchedConds.begin(); iter != rcRefSchedConds.end(); iter++)
		{			
			IloNumVar & rcThisVar = rcAbsDiffVars[rcAbsDiffVars.getSize() - 1];
			IloNumVar & rcThisSchedCondTypeVar = rcSchedCondTypeVars[iter->first];
			rcConst.add(IloRange(rcEnv,
				0.0,
				rcThisVar - (rcThisSchedCondTypeVar - iter->second),
				IloInfinity
				));
			rcConst.add(IloRange(rcEnv,
				0.0,
				rcThisVar - (iter->second - rcThisSchedCondTypeVar),
				IloInfinity
				));
		}
		IloExpr cObjExpr(rcEnv);	
		rcObj = IloMinimize(rcEnv, rcAbsDiffVars[0]);
		rcModel.add(rcObj);
		rcModel.add(rcConst);
	}

	virtual void GenRefGuidedSearchModel_MinSum(const IloEnv & rcEnv, SchedConds & rcRefSchedConds, SchedCondTypeVars & rcSchedCondTypeVars, IloModel & rcModel)
	{
		IloNumVarArray rcAbsDiffVars(rcEnv);
		IloRangeArray rcConst(rcEnv);
		IloObjective rcObj(rcEnv);
		rcAbsDiffVars = IloNumVarArray(rcEnv);
		//m_cUnschedCoreComputer.PrintSchedConds(rcRefSchedConds, cout); cout << endl; system("pause");
		for (typename SchedConds::iterator iter = rcRefSchedConds.begin(); iter != rcRefSchedConds.end(); iter++)
		{
			rcAbsDiffVars.add(IloNumVar(rcEnv, 0.0, IloInfinity));
			IloNumVar & rcThisVar = rcAbsDiffVars[rcAbsDiffVars.getSize() - 1];
			IloNumVar & rcThisSchedCondTypeVar = rcSchedCondTypeVars[iter->first];
			rcConst.add(IloRange(rcEnv,
				0.0,
				rcThisVar - (rcThisSchedCondTypeVar - iter->second),
				IloInfinity
			));
			rcConst.add(IloRange(rcEnv,
				0.0,
				rcThisVar - (iter->second - rcThisSchedCondTypeVar),
				IloInfinity
			));			
		}
		IloExpr cObjExpr(rcEnv);
		
		for (int i = 0; i < rcAbsDiffVars.getSize(); i++)
		{
			cObjExpr += rcAbsDiffVars[i];
		}		
		cObjExpr -= 100;
		rcObj = IloMinimize(rcEnv, cObjExpr);
		rcModel.add(rcObj);
		rcModel.add(rcConst);
	}

	virtual void GenReverseRefGuidedSearchModel(const IloEnv & rcEnv, SchedConds & rcRefSchedConds, SchedCondTypeVars & rcSchedCondTypeVars, IloModel & rcModel)
	{
		IloNumVarArray rcAbsDiffVars(rcEnv);
		IloRangeArray rcConst(rcEnv);
		IloObjective rcObj(rcEnv);
		rcAbsDiffVars = IloNumVarArray(rcEnv);
		rcAbsDiffVars.add(IloNumVar(rcEnv, 0.0, IloInfinity));
		for (typename SchedConds::iterator iter = rcRefSchedConds.begin(); iter != rcRefSchedConds.end(); iter++)
		{
			IloNumVar & rcThisVar = rcAbsDiffVars[rcAbsDiffVars.getSize() - 1];
			IloNumVar & rcThisSchedCondTypeVar = rcSchedCondTypeVars[iter->first];
			rcConst.add(IloRange(rcEnv,
				0.0,
				rcThisVar - (rcThisSchedCondTypeVar - iter->second),
				IloInfinity
				));
			rcConst.add(IloRange(rcEnv,
				0.0,
				rcThisVar - (iter->second - rcThisSchedCondTypeVar),
				IloInfinity
				));
		}
		IloExpr cObjExpr(rcEnv);
		rcObj = IloMinimize(rcEnv, rcAbsDiffVars[0]);
		rcModel.add(rcObj);
		rcModel.add(rcConst);
	}

	virtual void GenFixedSchedCondsConst(SchedConds & rcFixedSchedConds, SchedCondTypeVars & rcSchedCondTypeVars, IloRangeArray & rcConst)
	{
		const IloEnv & rcEnv = rcConst.getEnv();
		for (typename SchedConds::iterator iterSchedCond = rcFixedSchedConds.begin(); iterSchedCond != rcFixedSchedConds.end(); iterSchedCond++)
		{
			if (rcSchedCondTypeVars.count(iterSchedCond->first) == 0) continue;
			IloNumVar & rcThisVar = rcSchedCondTypeVars[iterSchedCond->first];
			if (m_cUnschedCoreComputer.SchedCondMonotonicity(iterSchedCond->first) == UnschedCoreComputer_MonoInt_T::Up)
			{
				rcConst.add(IloRange(rcEnv, 
					-IloInfinity, 
					rcThisVar - iterSchedCond->second, 
					0.0
					));
			}
			else
			{
				rcConst.add(IloRange(rcEnv,
					-IloInfinity,
					iterSchedCond->second - rcThisVar,
					0.0
					));
			}
		}
	}
	
	virtual void BestFeasibleHook(bool & rbCanExit)
	{
		//For implementing hierarchical schedulability analysis
		return;

	}
};

