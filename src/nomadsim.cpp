#include "myutils.h"
#include "nomadsim.h"

static uint GetRand(uint Lo, uint Hi)
	{
	uint Min = min(Lo, Hi);
	uint Max = max(Lo, Hi);
	uint Range = Max - Min + 1;
	uint r = randu32()%Range;
	uint n = Min + r;
	asserta(n >= Min && n <= Max);
	return n;
	}

uint NomadSim::GetDistTargetCount(uint DistIndex) const
	{
	asserta(DistIndex < SIZE(m_DistVec));
	uint n = SIZE(m_DistVec[DistIndex]);
	return n;
	}

void NomadSim::GetRandomNloNhi(uint &Nlo, uint &Nhi) const
	{
	uint r1 = GetRand(m_MinNlo, m_MaxNhi);
	uint r2 = GetRand(m_MinNlo, m_MaxNhi);
	Nlo = min(r1, r2);
	Nhi = max(r1, r2);
	}

uint NomadSim::GetRandomDistIndex() const
	{
	uint DistCount = SIZE(m_DistVec);
	asserta(DistCount > 0);
	uint DistIndex = randu32()%DistCount;
	return DistIndex;
	}

void NomadSim::GetRandomDistIndexes(uint DistCount,
  vector<uint> &DistIndexes) const
	{
	const uint Size = SIZE(m_DistVec);
	asserta(DistCount <= Size);
	DistIndexes.clear();
	for (uint i = 0; i < Size; ++i)
		DistIndexes.push_back(i);
	random_shuffle(DistIndexes.begin(), DistIndexes.end());
	DistIndexes.resize(DistCount);
	}

void NomadSim::GenerateCountsVec1(uint DistIndex, 
  uint SampleCount, vector<vector<uint> > &CountsVec) const
	{
	CountsVec.clear();

	uint Nlo, Nhi;
	GetRandomNloNhi(Nlo, Nhi);
	asserta(Nhi >= Nlo);

	if (DistIndex == UINT_MAX)
		DistIndex = GetRandomDistIndex();
	asserta(DistIndex < SIZE(m_DistVec));
	const vector<uint> &Dist = m_DistVec[DistIndex];

	vector<uint> DopeVec;
	const uint TargetCount = SIZE(Dist);
	CountsVec.resize(TargetCount);
	for (uint i = 0; i < TargetCount; ++i)
		{
		CountsVec[i].resize(SampleCount);

		uint n = Dist[i];
		for (uint j = 0; j < n; ++j)
			DopeVec.push_back(i);
		}
	const uint D = SIZE(DopeVec);

	for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		uint N = Nlo + randu32()%(Nhi - Nlo + 1);
		for (uint i = 0; i < N; ++i)
			{
			uint r = randu32()%D;
			uint TargetIndex = DopeVec[r];
			CountsVec[TargetIndex][SampleIndex] += 1;
			}
		}
	}

void NomadSim::GetRandomSampleCounts(uint DistCount, uint TotalSamples,
  vector<uint> &SampleCounts) const
	{
	SampleCounts.clear();

	uint FirstPct = 40 + randu32()%40;
	uint FirstCount = (TotalSamples*FirstPct)/100;
	if (FirstCount > TotalSamples - DistCount)
		FirstCount = TotalSamples - DistCount - 1;
	if (FirstCount == 0)
		FirstCount = 1;
	SampleCounts.push_back(FirstCount);

	uint RestCount = TotalSamples - FirstCount;
	for (uint i = 1; i + 1 < DistCount; ++i)
		{
		asserta(RestCount > 0);
		uint n = RestCount/DistCount;
		if (n == 0)
			n = 1;
		RestCount -= n;
		SampleCounts.push_back(n);
		}
	SampleCounts.push_back(RestCount);

	uint SumSamples2 = 0;
	for (uint i = 0; i < DistCount; ++i)
		SumSamples2 += SampleCounts[i];
	asserta(SumSamples2 == TotalSamples);

	uint T10 = TotalSamples/10;
	if (T10 == 0)
		T10 = 1;
	for (uint i = 0; i < 10; ++i)
		{
		uint d1 = randu32()%DistCount;
		if (d1 == 0)
			d1 = 1;
		uint d2 = randu32()%DistCount;
		if (d2 == 0)
			d2 = 1;
		uint n = randu32()%T10 + 1;
		if (SampleCounts[d1] > n )
			{
			SampleCounts[d1] -= n;
			SampleCounts[d2] += n;
			}
		}

	uint SumSamples3 = 0;
	for (uint i = 0; i < DistCount; ++i)
		SumSamples3 += SampleCounts[i];
	asserta(SumSamples3 == TotalSamples);
	}

void NomadSim::GenerateCountsVecN(uint DistCount,
  uint SampleCount, vector<vector<uint> > &CountsVec) const
	{
	CountsVec.clear();

	vector<uint> DistIndexes;
	GetRandomDistIndexes(DistCount, DistIndexes);

	uint TargetCount = 0;
	for (uint i = 0; i < DistCount; ++i)
		{
		uint DistIndex = DistIndexes[i];
		uint TargetCount_i = GetDistTargetCount(DistIndex);
		TargetCount = max(TargetCount, TargetCount_i);
		}
	asserta(TargetCount > 0);
	CountsVec.resize(TargetCount);

	vector<uint> SampleCounts;
	GetRandomSampleCounts(DistCount, SampleCount, SampleCounts);
	uint SumSampleCount = 0;
	for (uint i = 0; i < DistCount; ++i)
		{
		uint n = SampleCounts[i];
		SumSampleCount += n;
		}
	asserta(SumSampleCount == SampleCount);
	CountsVec.resize(TargetCount);
	for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		CountsVec[TargetIndex].resize(SampleCount);

	uint SumSampleCount2 = 0;
	uint SumTargetCount2 = 0;
	for (uint i = 0; i < DistCount; ++i)
		{
		uint DistIndex_i = DistIndexes[i];
		uint SampleCount_i = SampleCounts[i];

		uint TargetCount_i = GetDistTargetCount(DistIndex_i);

		vector<vector<uint> > CountsVec_i;
		GenerateCountsVec1(DistIndex_i, SampleCount_i, CountsVec_i);
		asserta(SIZE(CountsVec_i) == TargetCount_i);

		for (uint TargetIndex = 0; TargetIndex < TargetCount_i; ++TargetIndex)
			{
			for (uint SampleIndex = 0; SampleIndex < SampleCount_i; ++SampleIndex)
				{
				uint n = CountsVec_i[TargetIndex][SampleIndex];

				uint SampleIndex2 = SampleIndex + SumSampleCount2;
				asserta(TargetIndex < SIZE(CountsVec));
				asserta(SampleIndex2 < SIZE(CountsVec[TargetIndex]));
				CountsVec[TargetIndex][SampleIndex2] = n;
				}
			}
		SumSampleCount2 += SampleCount_i;
		}
	}
