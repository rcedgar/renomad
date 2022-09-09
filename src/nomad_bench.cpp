#include "myutils.h"
#include "params.h"
#include "mxreader.h"
#include "nomad.h"

static FILE *g_fOut;

void InitNomad(uint ka, uint kt, const vector<string> &SampleNames);

bool g_IsPos;
uint g_TP;
uint g_TN;
uint g_FP;
uint g_FN;

void InitBench()
	{
	g_TP = 0;
	g_TN = 0;
	g_FP = 0;
	g_FN = 0;
	}

void ReportBench(const Nomad &ND)
	{
	if (g_fOut == 0 && optset_output)
		g_fOut = CreateStdioFile(opt(output));

	double Recall = g_TP/double(g_TP + g_FN);
	double FDR = g_FP/double(g_FP + g_TP);

	Progress("\n");
	ProgressLog("Recall %.4f, FDR %.4f", Recall, FDR);
	ProgressLog(" MAT=%u", g_MinAnchorTotalBeforeFiltering);
	ProgressLog(" MSA=%u", g_MinSamplesPerAnchor);
	ProgressLog(" MTA=%u", g_MinTargetsPerAnchor);
	ProgressLog(" MST=%u", g_MinSampleTotal);
	ProgressLog(" MTMC=%u", g_MinTargetMaxCount);
	ProgressLog(" SC=%u", g_SplitCount);
	ProgressLog(" FC=%u", g_FuncCount);
	ProgressLog(" Try=%u", g_RetryCount);
	ProgressLog(" Bon=%c", tof(g_DoBonferroni));
	ProgressLog(" tp=%u", g_TP);
	ProgressLog(" tn=%u", g_TN);
	ProgressLog(" fp=%u", g_FP);
	ProgressLog(" fn=%u", g_FN);
	ProgressLog("\n");

	fprintf(g_fOut, "Recall=%.4f\tFDR=%.4f", Recall, FDR);
	fprintf(g_fOut, "\tMAT=%u", g_MinAnchorTotalBeforeFiltering);
	fprintf(g_fOut, "\tMSA=%u", g_MinSamplesPerAnchor);
	fprintf(g_fOut, "\tMTA=%u", g_MinTargetsPerAnchor);
	fprintf(g_fOut, "\tMST=%u", g_MinSampleTotal);
	fprintf(g_fOut, "\tMTMC=%u", g_MinTargetMaxCount);
	fprintf(g_fOut, "\tSC=%u", g_SplitCount);
	fprintf(g_fOut, "\tFC=%u", g_FuncCount);
	fprintf(g_fOut, "\tTry=%u", g_RetryCount);
	fprintf(g_fOut, "\tBon=%c", tof(g_DoBonferroni));
	fprintf(g_fOut, "\ttp=%u", g_TP);
	fprintf(g_fOut, "\ttn=%u", g_TN);
	fprintf(g_fOut, "\tfp=%u", g_FP);
	fprintf(g_fOut, "\tfn=%u", g_FN);
	fprintf(g_fOut, "\t\n");
	}

static bool OnStart(uint ka, uint kt,
  const vector<string> &SampleNames, uint64 ptrNd)
	{
	Nomad &Nd = *(Nomad *) ptrNd;
	Nd.Init(ka, kt, SampleNames, 1000);
	return true;
	}

static bool OnAnchor(const Atab &A, uint64 ptrNd)
	{
	Nomad &Nd = *(Nomad *) ptrNd;
	double P = Nd.CalcP();
	if (g_IsPos)
		{
		if (P <= g_MaxP)
			++g_TP;
		else
			++g_FN;
		}
	else
		{
		if (P <= g_MaxP)
			++g_FP;
		else
			++g_TN;
		}
	return true;
	}

static void BenchMx(const string &FileName, Nomad &Nd)
	{
	Progress("%s\n", FileName.c_str());
	MxReader MR;
	uint64 ptrNd = uint64(&Nd);
	MR.Run(FileName, Nd, 27, 27, OnStart, OnAnchor, ptrNd);
	}

static void BenchParams()
	{
	Nomad Nd;

	InitBench();

	g_IsPos = true;
	BenchMx("../simx/pos.sam10", Nd);
	BenchMx("../simx/pos.sam100", Nd);
	BenchMx("../simx/pos.sam1000", Nd);

	g_IsPos = false;
	BenchMx("../simx/neg.sam10", Nd);
	BenchMx("../simx/neg.sam100", Nd);
	BenchMx("../simx/neg.sam1000", Nd);

	ReportBench(Nd);
	}

void cmd_nomad_bench()
	{
	opt(nomad_bench);

	{
	for (uint Iter = 0; Iter < 4; ++Iter)
		{
		NomadParams NP1;
		SetParams_Struct(NP1);
		BenchParams();
		}
	}

	{
	NomadParams NP;
	NP.RetryCount = 2;
	SetParams_Struct(NP);
	BenchParams();
	}

	{
	NomadParams NP;
	NP.RetryCount = 4;
	SetParams_Struct(NP);
	BenchParams();
	}

	{
	NomadParams NP;
	NP.RetryCount = 8;
	SetParams_Struct(NP);
	BenchParams();
	}

	{
	NomadParams NP;
	NP.RetryCount = 16;
	SetParams_Struct(NP);
	BenchParams();
	}

	{
	NomadParams NP;
	NP.DoBonferroni = false;
	SetParams_Struct(NP);
	BenchParams();
	}

	{
	NomadParams NP;
	NP.MinSampleTotal = 8;
	SetParams_Struct(NP);
	BenchParams();
	}

	{
	NomadParams NP;
	NP.MinSampleTotal = 16;
	SetParams_Struct(NP);
	BenchParams();
	}

	{
	NomadParams NP;
	NP.MinTargetMaxCount = 2;
	SetParams_Struct(NP);
	BenchParams();
	}

	{
	NomadParams NP;
	NP.MinTargetMaxCount = 4;
	SetParams_Struct(NP);
	BenchParams();
	}

	{
	NomadParams NP;
	NP.MinTargetMaxCount = 8;
	SetParams_Struct(NP);
	BenchParams();
	}

	{
	NomadParams NP2;
	NP2.SplitCount = 100;
	SetParams_Struct(NP2);
	BenchParams();
	}

	{
	NomadParams NP3;
	NP3.SplitCount = 200;
	SetParams_Struct(NP3);
	BenchParams();
	}

	{
	NomadParams NP4;
	NP4.FuncCount = 20;
	SetParams_Struct(NP4);
	BenchParams();
	}

	{
	NomadParams NP5;
	NP5.FuncCount = 40;
	SetParams_Struct(NP5);
	BenchParams();
	}

	{
	NomadParams NP6;
	NP6.FuncCount = 5;
	SetParams_Struct(NP6);
	BenchParams();
	}
	}
