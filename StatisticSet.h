#pragma once
#include <iostream>
using namespace std;
#include <set>
#include <cstring>
#include <string.h>

#define ITEMNAME_LEN	256
#define ITEMSTRING_LEN	512
#define ITEMTYPE_VAL	0
#define ITEMTYPE_STR	1
class StatisticSet
{
public:
	class StatisticItem
	{	
	private:
		char m_axItemName[ITEMNAME_LEN];
		double m_dValue;
		char m_axString[ITEMSTRING_LEN];
		int m_iType;
	public:
		StatisticItem();
		StatisticItem(const char axName[], double dValue);
		StatisticItem(const char axName[], const char axString[]);
		~StatisticItem();		
		int getType() const	{ return m_iType; }
		double getValue() const	{ return m_dValue; }
		const char * getString() const	{ return m_axString; };
		const char * getItemName() const	{ return m_axItemName; }
		friend bool operator < (const StatisticItem & rcRhs, const StatisticItem & rcLhs);
		friend std::ostream & operator<<(std::ostream& os, const StatisticItem & obj);
	};
private:
	set<StatisticItem> m_setValueItems;//may be in the future we will have other type of items.
public:
	StatisticSet();
	~StatisticSet();	
	StatisticItem getItem(const char axName[]);
	void Clear();
	bool IsItemExist(char axName[]);
	bool setItem(const char axName[], double dValue);
	bool setItem(const char axName[], int iValue);
	void IncItem(const char axName[], double dIncValue);
	bool setItem(const char axName[], char axString[]);
	void WriteStatisticImage(char axFileName[]);
	void ReadStatisticImage(char axFileName[]);
	void WriteStatisticText(char axFileName[]);	
	void DisplayStatistic();
};

