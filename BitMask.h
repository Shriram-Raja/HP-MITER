#pragma once
#include <vector>
#include <fstream>
#include <unordered_map>
#include <map>
#include <set>
using namespace std;
//interesting tool class
class BitMask {
public:

protected:
	vector<unsigned int> m_vectorData;
	int m_iRightMostOne;
	int m_iLeftMostOne;
	int m_iNumOnes;
	int m_iWidth;	
public:
	BitMask():m_vectorData(1), m_iWidth(1) {};
	BitMask(int iWidth);
	BitMask(const char axStr[]);
	void set(int iPos);	
	int get(int iPos) const;
	inline int getWidth() const { return m_iWidth; }
	inline int getRightMostPos() const { return m_iRightMostOne; }
	inline int getLeftMostPos() const { return m_iLeftMostOne; }
	inline int getNumOnes() const { return m_iNumOnes; }
	bool operator >= (const BitMask & rhs) const;
	bool operator <= (const BitMask & rhs) const;
	bool operator < (const BitMask & rhs) const;
protected:
};

class BitMaskBasedSubSupsetQuery {
public:

protected:
	typedef map<int, set<BitMask>> RightMostOnePartition;
	typedef map<int, RightMostOnePartition> NumOnesPartition;
	NumOnesPartition m_cPartition;	
public:
	BitMaskBasedSubSupsetQuery() {}
	void insert(const BitMask & rcBitMask);
	void insert(BitMaskBasedSubSupsetQuery & rcOther);
	pair<bool, BitMask> existSupset(BitMask & rcBitMask);
	pair<bool, BitMask> existSubset(BitMask & rcBitMask);
protected:
};

void testBitMask();