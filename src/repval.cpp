#include "myutils.h"
#include "sort.h"

static uint g_DopeVecScale = 32;
static uint g_REPIters = 100;
static uint g_TN = 5;

#define TRACE	1

void LogCountsVec(const vector<vector<uint> > &CountsVec)
	{
	const uint TargetCount = SIZE(CountsVec);
	if (TargetCount == 0)
		{
		Log("Empty CountsVec\n");
		return;
		}
	const uint SampleCount = SIZE(CountsVec[0]);
	Log("%8.8s", "Tgt");
	for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		string s;
		Ps(s, "Sam%u", SampleIndex);
		Log(" %5.5s", s.c_str());
		}
	Log("\n");

	for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		string s;
		Ps(s, "T%u", TargetIndex);
		Log("%8.8s", s.c_str());
		for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			uint Count = CountsVec[TargetIndex][SampleIndex];
			Log(" %5u", Count);
			}
		Log("\n");
		}
	}

void LogFreqMx(const vector<vector<double> > &FreqMx)
	{
	const uint TargetCount = SIZE(FreqMx);
	if (TargetCount == 0)
		{
		Log("Empty FreqMx\n");
		return;
		}
	const uint SampleCount = SIZE(FreqMx[0]);
	Log("%8.8s", "Tgt");
	for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		string s;
		Ps(s, "Sam%u", SampleIndex);
		Log("  %8.8s", s.c_str());
		}
	Log("\n");

	for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		string s;
		Ps(s, "T%u", TargetIndex);
		Log("%8.8s", s.c_str());
		for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			double Freq = FreqMx[TargetIndex][SampleIndex];
			Log("  %8.6f", Freq);
			}
		Log("\n");
		}
	}

static double CalcSumFreqDiff(const vector<vector<double> > &TargetFreqMx,
  uint SampleIndex1, uint SampleIndex2)
	{
#if 0//TRACE
	Log("\n");
	Log("CalcSumFreqDiff(%u, %u)\n", SampleIndex1, SampleIndex2);
#endif

	const uint TargetCount = SIZE(TargetFreqMx);
	double SumFreqDiff = 0;
	for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		double Freq1 = TargetFreqMx[TargetIndex][SampleIndex1];
		double Freq2 = TargetFreqMx[TargetIndex][SampleIndex2];
		double Diff = fabs(Freq1 - Freq2);
		SumFreqDiff += Diff;
#if 0//TRACE
		Log(" T%u %.3g-%.3g", TargetIndex, Freq1, Freq2);
#endif
		}
#if 0//TRACE
	Log(", Sum=%.3g\n", SumFreqDiff);
#endif
	return SumFreqDiff;
	}

static double CalcTestStat(const vector<vector<double> > &TargetFreqMx)
	{
	const uint TargetCount = SIZE(TargetFreqMx);
	asserta(TargetCount > 0);
	const uint SampleCount = SIZE(TargetFreqMx[0]);
	vector<double> FDs;
	for (uint SampleIndex1 = 0; SampleIndex1 < SampleCount; ++SampleIndex1)
		{
		for (uint SampleIndex2 = SampleIndex1 + 1; 
		  SampleIndex2 < SampleCount; ++SampleIndex2)
			{
			double FD = CalcSumFreqDiff(TargetFreqMx, SampleIndex1, SampleIndex2);
			FDs.push_back(FD);
			}
		}
	const uint PairCount = SIZE(FDs);

	QuickSortInPlaceDesc(FDs.data(), PairCount);

	double T1 = double(PairCount)/g_TN;
	double T2 = ceil(T1 + 0.5);
	uint T = uint(T2);
	if (T == 0)
		T = 1;

//	double Stat = FDs[T-1];
	double Stat = FDs[0];

#if TRACE
	Log("\n");
	Log("CalcTestStat()\n");
	Log("FDs ");
	for (uint PairIndex = 0; PairIndex < SIZE(FDs); ++PairIndex)
		{
		if (PairIndex > 0 && PairIndex%20 == 0)
			Log("\n");
		Log("  %.3g", FDs[PairIndex]);
		}
	Log("\nStat=%.3g\n", Stat);
	Log("\n");
#endif
	return Stat;
	}

static void GetSimFreqMx(uint TargetCount, const vector<uint> &SampleTotals,
  const vector<uint> &DopeVec, vector<vector<double> > &SimFreqMx)
	{
	const uint SampleCount = SIZE(SampleTotals);
	SimFreqMx.clear();
	SimFreqMx.resize(TargetCount);
	for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		SimFreqMx[TargetIndex].resize(SampleCount);

	uint D = SIZE(DopeVec);
	asserta(D > 0);

	for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		uint SampleTotal = SampleTotals[SampleIndex];
		if (SampleTotal == 0)
			continue;
		double Sliver = 1.0/SampleTotal;
		for (uint i = 0; i < SampleTotal; ++i)
			{
			uint r = randu32()%D;
			uint TargetIndex = DopeVec[r];
			asserta(TargetIndex < TargetCount);
			SimFreqMx[TargetIndex][SampleIndex] += Sliver;
			}
		}
	}

double CalcREPvalue(const vector<vector<uint> > &CountsVec)
	{
	const uint TargetCount = SIZE(CountsVec);
	if (TargetCount <= 1)
		return 1.0;

	const uint SampleCount = SIZE(CountsVec[0]);

#if TRACE
	Log("\n");
	Log("CalcREPvalue()\n");
	LogCountsVec(CountsVec);
#endif

	vector<uint> TargetTotals(TargetCount);
	vector<uint> SampleTotals(SampleCount);
	uint AnchorTotal = 0;
	for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			uint Count = CountsVec[TargetIndex][SampleIndex];
			AnchorTotal += Count;
			TargetTotals[TargetIndex] += Count;
			SampleTotals[SampleIndex] += Count;
			}
		}

	if (AnchorTotal == 0)
		return 1.0;

	double SumFreqs = 0;
	vector<uint> DopeVec;
	vector<double> TargetFreqs;
	for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		uint TargetTotal = TargetTotals[TargetIndex];
		double TargetFreq = double(TargetTotal)/double(AnchorTotal);
		SumFreqs += TargetFreq;
		TargetFreqs.push_back(TargetFreq);

		uint n = uint(TargetFreq*g_DopeVecScale + 0.5);
		for (uint i = 0; i < n; ++i)
			DopeVec.push_back(TargetIndex);
		}
	asserta(feq(SumFreqs, 1.0));
	if (DopeVec.empty())
		return 1.0;

	vector<vector<double> > TargetFreqMx(TargetCount);
	for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		TargetFreqMx[TargetIndex].resize(SampleCount);
		for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			const uint SampleTotal = SampleTotals[SampleIndex];
			uint n = CountsVec[TargetIndex][SampleIndex];
			asserta(n <= SampleTotal);
			double Freq = double(n)/double(SampleTotal);
			TargetFreqMx[TargetIndex][SampleIndex] = Freq;
			}
		}

	double TestStat = CalcTestStat(TargetFreqMx);
#if TRACE
	Log("\n");
	Log("FreqMx, TestStat %.3g\n", TestStat);
	LogFreqMx(TargetFreqMx);
#endif

	vector<vector<double> > SimFreqMx(TargetCount);
	for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		SimFreqMx[TargetIndex].resize(SampleCount);

	uint HitCount = 0;
	for (uint Iter = 0; Iter < g_REPIters; ++Iter)
		{
		GetSimFreqMx(TargetCount, SampleTotals, DopeVec, SimFreqMx);
		double SimTestStat = CalcTestStat(SimFreqMx);
		if (SimTestStat >= TestStat)
			++HitCount;
#if TRACE
		Log("\n");
		char op = '?';
		if (SimTestStat > TestStat)
			op = '>';
		else if (SimTestStat < TestStat)
			op = '<';
		else
			{
			asserta(SimTestStat == TestStat);
			op = '=';
			}
		Log("Iter %u, SimFreqMx, SimTestStat %.3g %c %.3g\n",
		  Iter, SimTestStat, op, TestStat);
		LogFreqMx(SimFreqMx);
#endif
		}
	double P = 1.0 - double(HitCount)/g_REPIters;
#if TRACE
	Log("P = %.3g\n", P);
#endif
	return P;
	}
