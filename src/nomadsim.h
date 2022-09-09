#pragma once

class NomadSim
	{
public:
	uint m_MinNlo = 1;
	uint m_MaxNhi = 100;

	uint m_MinSplitPct = 20;

	vector<vector<uint> > m_DistVec;

public:
	uint GetDistTargetCount(uint DistIndex) const;
	void GetRandomNloNhi(uint &Nlo, uint &Nhi) const;
	uint GetRandomDistIndex() const;
	void GetRandomDistIndexes(uint N, 
	  vector<uint> &DistIndexes) const;
	void GetRandomSampleCounts(uint DistCount, uint TotalSamples,
	  vector<uint> &SampleCounts) const;

	void GenerateCountsVec1(uint DistIndex, uint SampleCount,
	  vector<vector<uint> > &CountsVec) const;

	void GenerateCountsVecN(uint DistCount, uint SampleCount,
	  vector<vector<uint> > &CountsVec) const;
	};
