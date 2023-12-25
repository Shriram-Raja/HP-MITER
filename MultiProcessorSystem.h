#pragma once
#include <iostream>
using namespace std;
#include <vector>
namespace MPScheduling {
	//typedef _int64 ValueType;
	typedef int ValueType;
	class Task {
	public:

	protected:
		ValueType m_dWCET;
		ValueType m_dPeriod;
		ValueType m_dDeadline;
	public:
		Task();
		Task(ValueType dWCET, ValueType dPeriod, ValueType dDeadline);
		~Task() {}
		inline ValueType C() const { return m_dWCET; }
		inline ValueType T() const { return m_dPeriod; }
		inline ValueType D() const { return m_dDeadline; }
		double U() const { return (double)m_dWCET / (double)m_dPeriod; }
	protected:		
	};

	class MultiProcessorSystem {
	public:

	protected:
		vector<Task> m_vectorTasks;
		int m_iNumProcessors;
	public:
		MultiProcessorSystem() : m_iNumProcessors(0) {}
		MultiProcessorSystem(int iNumProcessors);
		~MultiProcessorSystem();
		const vector<Task> getTasks() { return m_vectorTasks; }
		inline const Task & getTask(int iIndex) { return m_vectorTasks.at(iIndex); }
		inline int getNumProcessors() const { return m_iNumProcessors; }
		inline int getNumTasks() const { return m_vectorTasks.size(); }
		void AddTask(Task rcTask) { m_vectorTasks.push_back(rcTask); }
		void Write(const char axFileName[]);
		void Read(const char axFileName[]);
		double TotalUtil() const;
	protected:
		void Write(ofstream & cOutputFile);
		void Read(ifstream & cInputFile);
	};

	class PriorityAssignment {
	public:

	protected:
		vector<int> m_vectorTask2Priority;
		vector<int> m_vectorPriority2Task;
	public:
		PriorityAssignment() {}
		PriorityAssignment(MultiProcessorSystem & rcSystem);
		int getPriority(int iTaskIndex) const;
		int getTask(int iPriority) const;
		int getPriorityByTask(int iTaskIndex) const { return getPriority(iTaskIndex); }
		int getTaskByPriority(int iPriority) const { return getTask(iPriority); }
		void CopyFrom_Strict(PriorityAssignment & rcOther) { *this = rcOther; }
		int getSize() const;
		void setPriority(int iTaskIndex, int iPriority);
		void unset(int iTaskIndex);
		void Write(const char axFileName[]);
		void Read(const char axFileName[]);
		void GenerateDkCPA(MultiProcessorSystem & rcSystem);
		void GenerateDMPA(MultiProcessorSystem & rcSystem);
		void GenerateDCMPA(MultiProcessorSystem & rcSystem);
		int GenerateHP_MUTER(MultiProcessorSystem & rcSystem, float THRESH, bool DkC_if_not_MUTER, bool Lateness);
	protected:
	};

	class RandomSystemGenerator {
	public:

	protected:

	public:
		RandomSystemGenerator();
		~RandomSystemGenerator() {}
		MultiProcessorSystem GenerateSystem(int iTaskNum, int iProcessorNum, double dTotalUtil, ValueType cPeriodLB, ValueType cPeriodUB);
		MultiProcessorSystem GenerateSystemConstrainedDeadline(int iTaskNum, int iProcessorNum, double dTotalUtil, ValueType cPeriodLB, ValueType cPeriodUB);
		MultiProcessorSystem GenerateSystemMLforRT(int iTaskNum, int iProcessorNum, double dTotalUtil, ValueType cPeriodLB, ValueType cPeriodUB);
	protected:
		vector<double> UUnifastDiscard(int iTaskNum, double dTotalUtil);
		vector<double> UUnifastMLforRT(int iTaskNum, double dTotalUtil);
		vector<ValueType> GeneratePeriod(int iTaskNum, ValueType lb, ValueType ub);
		vector<ValueType> GeneratePeriodLogUniform(int iTaskNum, ValueType lb, ValueType ub);		
	};
}
