#include "myutils.h"
#include "nomadsim.h"
#include "nomad.h"
#include "params.h"

void InitBench();
void ReportBench(const Nomad &Nd);

static uint g_Iters;

extern uint g_TP;
extern uint g_FP;
extern uint g_TN;
extern uint g_FN;
extern bool g_IsPos;

static NomadSim g_NS;
static Nomad g_Nomad;
static uint g_CurrSampleCount;

static void AddDist(NomadSim &g_NS, const string &s)
	{
	vector<string> Fields;
	Split(s, Fields, ',');
	vector<uint> Dist;
	const uint N = SIZE(Fields);
	for (uint i = 0; i < N; ++i)
		{
		uint n = StrToUint(Fields[i]);
		Dist.push_back(n);
		}
	g_NS.m_DistVec.push_back(Dist);
	}

static void AddDists(NomadSim &g_NS)
	{
	AddDist(g_NS, "100");
	AddDist(g_NS, "50,50");
	AddDist(g_NS, "25,25,50");
	AddDist(g_NS, "25,50,25");
	AddDist(g_NS, "10,20,70");
	AddDist(g_NS, "33,33,34");
	AddDist(g_NS, "10,10,10,80");
	AddDist(g_NS, "80,5,5,5,5");
	}

static void Do1(uint SampleCount,
  const vector<vector<uint> > &CountsVec,
  bool IsPos)
	{
	if (SampleCount != g_CurrSampleCount)
		{
		g_Nomad.Init(27, 27, SampleCount, 1000);
		g_CurrSampleCount = SampleCount;
		}
	g_Nomad.FromVec(CountsVec);
	double P = g_Nomad.CalcP();
	bool Hit = (P <= g_MaxP);
	if (IsPos)
		{
		if (Hit)
			++g_TP;
		else
			++g_FN;
		}
	else
		{
		if (Hit)
			++g_FP;
		else
			++g_TN;
		}
//	Log("Pos=%c P=%.3g\n", tof(IsPos), P);
	}

static void DoNeg(uint SampleCount)
	{
	vector<vector<uint> > CountsVec;
	uint DistIndex = g_NS.GetRandomDistIndex();
	g_NS.GenerateCountsVec1(DistIndex, SampleCount, CountsVec);
	Do1(SampleCount, CountsVec, false);
	}

static void DoPos(uint SampleCount)
	{
	vector<vector<uint> > CountsVec;
	// uint DistCount = 2 + randu32()%5;
	uint DistCount = 2;
	g_NS.GenerateCountsVecN(DistCount, SampleCount, CountsVec);
	Do1(SampleCount, CountsVec, true);
	}

static void Bench1(const NomadParams &NP)
	{
	SetParams_Struct(NP);
	InitBench();

	uint SampleCount = 100;
	for (uint Iter = 0; Iter < g_Iters; ++Iter)
		{
		if (Iter == g_Iters/3)
			SampleCount = 30;
		if (Iter == 2*g_Iters/3)
			SampleCount = 10;
		ProgressStep(Iter, g_Iters, "Working (%u)", SampleCount);
		if (Iter%10 == 0)
			DoPos(SampleCount);
		else
			DoNeg(SampleCount);
		}

	ReportBench(g_Nomad);
	}
	
void cmd_nomad_sim_bench()
	{
	g_Iters = StrToUint(opt(nomad_sim_bench));

	AddDists(g_NS);

	{
	for (uint Iter = 0; Iter < 4; ++Iter)
		{
		NomadParams NP1;
		Bench1(NP1);
		}
	}

	for (uint Iter = 0; Iter < 4; ++Iter)
		{
		NomadParams NP;
		NP.RetryCount = 2;
		Bench1(NP);
		}

	for (uint Iter = 0; Iter < 4; ++Iter)
		{
		NomadParams NP;
		NP.RetryCount = 4;
		Bench1(NP);
		}

	for (uint Iter = 0; Iter < 4; ++Iter)
		{
		NomadParams NP;
		NP.RetryCount = 8;
		Bench1(NP);
		}

	for (uint Iter = 0; Iter < 4; ++Iter)
		{
		NomadParams NP;
		NP.RetryCount = 16;
		Bench1(NP);
		}

	{
	NomadParams NP;
	NP.DoBonferroni = false;
	Bench1(NP);
	}

	{
	NomadParams NP;
	NP.MinSampleTotal = 8;
	Bench1(NP);
	}

	{
	NomadParams NP;
	NP.MinSampleTotal = 16;
	Bench1(NP);
	}

	{
	NomadParams NP;
	NP.MinTargetMaxCount = 2;
	Bench1(NP);
	}

	{
	NomadParams NP;
	NP.MinTargetMaxCount = 4;
	Bench1(NP);
	}

	{
	NomadParams NP;
	NP.MinTargetMaxCount = 8;
	Bench1(NP);
	}

	{
	NomadParams NP2;
	NP2.SplitCount = 100;
	Bench1(NP2);
	}

	{
	NomadParams NP3;
	NP3.SplitCount = 200;
	Bench1(NP3);
	}

	{
	NomadParams NP4;
	NP4.FuncCount = 20;
	Bench1(NP4);
	}

	{
	NomadParams NP5;
	NP5.FuncCount = 40;
	Bench1(NP5);
	}

	{
	NomadParams NP6;
	NP6.FuncCount = 5;
	Bench1(NP6);
	}
	}
