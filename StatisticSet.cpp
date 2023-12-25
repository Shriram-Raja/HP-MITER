#include "stdafx.h"
#include "StatisticSet.h"
#include <algorithm>
#include <cmath>
#include "Util.h"
#include <fstream>
StatisticSet::StatisticItem::StatisticItem()
{
	memset(m_axItemName, 0, ITEMNAME_LEN);
}

StatisticSet::StatisticItem::StatisticItem(const char axName[], double dValue)
{
	my_assert(strlen(axName) <= ITEMNAME_LEN, "Name of the statistic item exceed max length");
	sprintf(m_axItemName, "%s", axName);
	m_dValue = dValue;
	m_iType = ITEMTYPE_VAL;
}

StatisticSet::StatisticItem::StatisticItem(const char axName[], const char axString[])
{
	my_assert(strlen(axName) <= ITEMNAME_LEN, "Name of the statistic item exceed max length");
	my_assert(strlen(axString) <= ITEMSTRING_LEN, "Length of the statistic item exceed max length");
	sprintf(m_axItemName, "%s", axName);
	sprintf(m_axString, "%s", axString);
	m_iType = ITEMTYPE_STR;
}

bool operator < (const StatisticSet::StatisticItem & rcRhs, const StatisticSet::StatisticItem & rcLhs)
{
	int iRhsLen = strlen(rcRhs.m_axItemName);
	int iLhsLen = strlen(rcLhs.m_axItemName);
	int iLen = min(iRhsLen, iLhsLen);
	for (int i = 0; i < iLen; i++)
	{
		if (rcRhs.m_axItemName[i] < rcLhs.m_axItemName[i])
		{
			return true;
		}
		else if (rcRhs.m_axItemName[i] > rcLhs.m_axItemName[i])
		{
			return false;
		}
	}
	return iRhsLen < iLhsLen;
}

std::ostream & operator<<(std::ostream& os, const StatisticSet::StatisticItem & obj)
{	
	if (obj.getType() == ITEMTYPE_STR)
		os << obj.m_axString;
	else if (obj.getType() == ITEMTYPE_VAL)
		os << obj.getValue();
	else
		my_assert(false, "Illegal Item Type Confronted");
	return os;
}

StatisticSet::StatisticItem::~StatisticItem()
{

}

StatisticSet::StatisticSet()
{
}


StatisticSet::~StatisticSet()
{

}

void StatisticSet::Clear()
{
	m_setValueItems.clear();
}

bool StatisticSet::IsItemExist(char axName[])
{
	set<StatisticItem>::iterator iterFind = m_setValueItems.find(StatisticItem(axName, (double)0));
	if (iterFind == m_setValueItems.end())
	{
		return false;
	}
	return true;
}

bool StatisticSet::setItem(const char axName[], double dValue)
{
	std::pair<set<StatisticItem>::iterator, bool> cPair = m_setValueItems.insert(StatisticItem(axName, dValue));
	if (!cPair.second)
	{
		m_setValueItems.erase(cPair.first);
		m_setValueItems.insert(StatisticItem(axName, dValue));
	}
	return cPair.second;
}

bool StatisticSet::setItem(const char axName[], int iValue)
{
	return setItem(axName, (double)iValue);
}

bool StatisticSet::setItem(const char axName[], char axString[])
{
	std::pair<set<StatisticItem>::iterator, bool> cPair = m_setValueItems.insert(StatisticItem(axName, axString));
	if (!cPair.second)
	{
		m_setValueItems.erase(cPair.first);
		m_setValueItems.insert(StatisticItem(axName, axString));
	}
	return cPair.second;
}

void StatisticSet::IncItem(const char axName[], double dIncValue)
{
	StatisticItem cItem = getItem(axName);
	my_assert(cItem.getType() == ITEMTYPE_VAL, "Algorithmic Error. Attempting to increse non-value item");
	setItem(axName, cItem.getValue() + dIncValue);

}

StatisticSet::StatisticItem StatisticSet::getItem(const char axName[])
{
	set<StatisticItem>::iterator iterFind = m_setValueItems.find(StatisticItem(axName, (double)0));
	if (iterFind == m_setValueItems.end())
	{
		cout << "Cannot Find " << axName << endl;
		my_assert(false, "Attempting to Find Nonexistent Statistic Item ");
	}
	return *iterFind;
}

void StatisticSet::WriteStatisticImage(char axFileName[])
{
	ofstream ofstreamOutputFile(axFileName, ios::out | ios::binary);
	int iValueItemSize = m_setValueItems.size();
	ofstreamOutputFile.write((const char *)&iValueItemSize, sizeof(int));
	for (set<StatisticItem>::iterator iter = m_setValueItems.begin();
		iter != m_setValueItems.end(); iter++)
	{
		ofstreamOutputFile.write(iter->getItemName(), ITEMNAME_LEN);
		int iType = iter->getType();
		ofstreamOutputFile.write((const char *)&iType, sizeof(int));
		if (iter->getType() == ITEMTYPE_VAL)
		{
			double dValue = iter->getValue();
			ofstreamOutputFile.write((const char *)&dValue, sizeof(double));
		}			
		else if (iter->getType() == ITEMTYPE_STR)
			ofstreamOutputFile.write((const char *)iter->getString(), ITEMSTRING_LEN);
		else
			my_assert(false, "Illegal Item Type Confronted");
	}
	ofstreamOutputFile.close();
	return;
}

void StatisticSet::ReadStatisticImage(char axFileName[])
{
	ifstream ifstreamInputFile(axFileName, ios::in | ios::binary);
	if (!ifstreamInputFile.is_open()) cout << axFileName << endl;
	my_assert(ifstreamInputFile.is_open(), "Cannot Find Statistic Item File");
	Clear();
	int iValueItemSize = 0;
	ifstreamInputFile.read((char *)&iValueItemSize, sizeof(int));
	char axItemName[ITEMNAME_LEN] = { 0 };
	char axItemString[ITEMSTRING_LEN] = { 0 };
	double dValue = 0;
	int iType = 0;
	m_setValueItems.clear();
	for (int i = 0; i < iValueItemSize; i++)
	{
		ifstreamInputFile.read(axItemName, ITEMNAME_LEN);
		ifstreamInputFile.read((char *)&iType, sizeof(int));
		if (iType == ITEMTYPE_VAL)
		{
			ifstreamInputFile.read((char *)&dValue, sizeof(double));
			m_setValueItems.insert(StatisticItem(axItemName, dValue));
		}
		else if (iType == ITEMTYPE_STR)
		{
			ifstreamInputFile.read(axItemString, ITEMSTRING_LEN);
			m_setValueItems.insert(StatisticItem(axItemName, axItemString));
		}
		else
			my_assert(false, "Illegal Item Type Confronted");
		
	}
	ifstreamInputFile.close();
}

void StatisticSet::WriteStatisticText(char axFileName[])
{
	ofstream ofstreamOutputFile(axFileName, ios::out);
	int iValueItemSize = m_setValueItems.size();
	ofstreamOutputFile << "Total " << iValueItemSize << " Items" << endl << endl;
	for (set<StatisticItem>::iterator iter = m_setValueItems.begin();
		iter != m_setValueItems.end(); iter++)
	{
		ofstreamOutputFile << "--------------------------------------" << endl;
		ofstreamOutputFile << iter->getItemName() << ": " ;
		if (iter->getType() == ITEMTYPE_VAL)
			ofstreamOutputFile << iter->getValue() << endl;
		else if (iter->getType() == ITEMTYPE_STR)
			ofstreamOutputFile << iter->getString() << endl;
		else
			my_assert(false, "Illegal Item Type Confronted");
	}
	ofstreamOutputFile.close();
	return;
}

void StatisticSet::DisplayStatistic()
{
	for (set<StatisticItem>::iterator iter = m_setValueItems.begin();
		iter != m_setValueItems.end(); iter++)
	{
		cout << "----------------------------" << endl;
		cout << iter->getItemName() << ":\t";
		if (iter->getType() == ITEMTYPE_VAL)
			cout << iter->getValue() << endl;
		else if (iter->getType() == ITEMTYPE_STR)
			cout << iter->getString() << endl;
		else
			my_assert(false, "Illegal Item Type Confronted");
	}
}