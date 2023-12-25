#include "stdafx.h"
#include "BitMask.h"
#include <cmath>
#include <assert.h>
#include <algorithm>
#include <assert.h>
#include "Util.h"
#include <algorithm>
#include <cstring>

const unsigned int aiBitMaskSingleBitMask[] = {
	0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80, 0x100, 0x200, 0x400, 0x800, 0x1000, 0x2000, 
	0x4000, 0x8000, 0x10000, 0x20000, 0x40000, 0x80000, 0x100000, 0x200000, 0x400000, 0x800000, 
	0x1000000, 0x2000000, 0x4000000, 0x8000000, 0x10000000, 0x20000000, 0x40000000, 0x80000000,
};


BitMask::BitMask(int iWidth):m_iWidth(iWidth), m_iRightMostOne(-1), m_iLeftMostOne(m_iWidth), m_iNumOnes(0){
	m_vectorData = vector<unsigned int>(iWidth / sizeof(unsigned int) / 8 + 1, 0);	
}

BitMask::BitMask(const char axStr[]) {
	int iLen = strlen(axStr);
	*this = BitMask(iLen);
	for(int i = 0; i < iLen; i++){
		assert(axStr[i] == '1' || axStr[i] == '0');
		if (axStr[i] == '1') set(i);				
	}
}

void BitMask::set(int iPos) {
	assert(iPos < m_iWidth);
	int iCellIndex = iPos / (sizeof(unsigned int) << 3);
	int iOffset = iPos % (sizeof(unsigned int) << 3);
	m_vectorData[iCellIndex] = m_vectorData[iCellIndex] | (1 << iOffset);
	m_iRightMostOne = max(m_iRightMostOne, iPos);
	m_iLeftMostOne = min(m_iLeftMostOne, iPos);
	m_iNumOnes++;
}

int BitMask::get(int iPos) const {
	assert(iPos < m_iWidth);
	int iCellIndex = iPos / (sizeof(unsigned int) << 3);
	int iOffset = iPos % (sizeof(unsigned int) << 3);
	return m_vectorData[iCellIndex] & ((1 << iOffset));
}

bool BitMask::operator >= (const BitMask & rhs) const {
	assert(rhs.getWidth() == getWidth());
	for (int i = 0; i < m_vectorData.size(); i++) 
		if ((m_vectorData[i] & rhs.m_vectorData[i]) != rhs.m_vectorData[i]) return false;
	return true;
}

bool BitMask::operator <= (const BitMask & rhs) const {
	assert(rhs.getWidth() == getWidth());
	for (int i = 0; i < m_vectorData.size(); i++)
		if ((m_vectorData[i] & rhs.m_vectorData[i]) != m_vectorData[i]) return false;
	return true;
}

bool BitMask::operator < (const BitMask & rhs) const {
	return less<vector<unsigned int>>()(this->m_vectorData, rhs.m_vectorData);
}

void BitMaskBasedSubSupsetQuery::insert(const BitMask & rcBitMask) {
	m_cPartition[rcBitMask.getNumOnes()][rcBitMask.getRightMostPos()].insert(rcBitMask);
}

void BitMaskBasedSubSupsetQuery::insert(BitMaskBasedSubSupsetQuery & rcOther) {
	for (auto & eleA : rcOther.m_cPartition) {
		for (auto & eleB : eleA.second) {
			for (auto & eleC : eleB.second) insert(eleC);
		}
	}
}

pair<bool, BitMask> BitMaskBasedSubSupsetQuery::existSupset(BitMask & rcBitMask) {
	if (m_cPartition[rcBitMask.getNumOnes()][rcBitMask.getRightMostPos()].count(rcBitMask)) return { true, rcBitMask };
	for (auto iterNumOnes = m_cPartition.lower_bound(rcBitMask.getNumOnes()); iterNumOnes != m_cPartition.end(); iterNumOnes++) {
		for (auto iterRightMost = iterNumOnes->second.lower_bound(rcBitMask.getRightMostPos()); iterRightMost != iterNumOnes->second.end(); iterRightMost++) {
			auto & rcCandidates = iterRightMost->second;
			for (auto & ele : rcCandidates) {
				if (rcBitMask <= ele) return { true, ele };
			}
		}
	}
	return {false, BitMask()};
}

pair<bool, BitMask> BitMaskBasedSubSupsetQuery::existSubset(BitMask & rcBitMask) {
	if (m_cPartition[rcBitMask.getNumOnes()][rcBitMask.getRightMostPos()].count(rcBitMask)) return { true, rcBitMask };
	for (auto iterNumOnes = m_cPartition.begin(); iterNumOnes != m_cPartition.end() && iterNumOnes->first <= rcBitMask.getNumOnes(); iterNumOnes++) {
		for (auto iterRightMost = iterNumOnes->second.begin(); iterRightMost != iterNumOnes->second.end() && iterRightMost->first <= rcBitMask.getRightMostPos(); iterRightMost++) {
			auto & rcCandidates = iterRightMost->second;
			for (auto & ele : rcCandidates) {
				if (rcBitMask >= ele) return { true, ele };
			}
		}		
	}
	return { false, BitMask() };
}

void testBitMask() {
	BitMask A("10110101"), B("10110111");	
	BitMaskBasedSubSupsetQuery cStore;
	cStore.insert(B);
	cout << cStore.existSubset(A).first << endl;
	cout << cStore.existSupset(A).first << endl;
}

