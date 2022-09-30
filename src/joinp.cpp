#include "myutils.h"
#include "joinera.h"
#include "nomad.h"
#include "params.h"

void ReadStringsFromFile(const string &FileName,
  vector<string> &Strs, bool DeleteEmptyStrings);
double CalcPvalueBrute(const vector<vector<uint> > &CountsVec);
void SetParams_CmdLine();

static uint g_HitCount = 0;
static FILE *g_fOut = 0;
static FILE *g_fTsv3Out = 0;
static uint g_LastMilDone = UINT_MAX;
static JoinerA *g_JA;

static bool OnAnchor(Atab &A, uint64 _notused_)
	{
	Nomad &Nd = (Nomad &) A;

	if (g_LastMilDone == UINT_MAX)
		{
		ProgressStep(0, 1001, "Scanning");
		g_LastMilDone = 0;
		}

	double P = Nd.CalcP();
	if (P <= g_MaxP)
		{
		double Pct = g_JA->GetPctDone();
		uint Mil = uint(Pct*10);
		ProgressStep(Mil, 1001, "%u significant anchors", ++g_HitCount);
		
		uint NS = Nd.m_SelectedSampleCount;
		uint NT = Nd.m_SelectedTargetCount;
		Pf(g_fOut, "%.2g\t%s\t%u\t%u\n",
		  P, A.m_Anchor, NS, NT);
		g_JA->m_A->ToFileTsv3(g_fTsv3Out);
		}
	return true;
	}

void cmd_joinp()
	{
	const string &InputFileName = opt(joinp);
	const string &OutputFileName = opt(output);
	const string &Tsv3FileName = opt(tsv3out);

	SetParams_CmdLine();
	
	g_fOut = CreateStdioFile(OutputFileName);
	g_fTsv3Out = CreateStdioFile(Tsv3FileName);

	vector<string> KDBNames;
	ReadStringsFromFile(InputFileName, KDBNames, true);
	const uint N = SIZE(KDBNames);

	vector<string> SampleNames;
	for (uint i = 0; i < N; ++i)
		{
		const string &KDBName = KDBNames[i];
		const char *p = strrchr(KDBName.c_str(), '/');
		if (p == 0)
			SampleNames.push_back(KDBName);
		else
			{
			string Name = string(p+1);
			SampleNames.push_back(Name);
			}
		}

	Nomad Nd;
	const uint ka = 27;
	const uint kt = 27;
	Nd.Init(ka, kt, SampleNames, 1000);

	g_HitCount = 0;
	g_LastMilDone = UINT_MAX;

	JoinerA JA;
	g_JA = &JA;
	JA.RunA(KDBNames, Nd, OnAnchor, 0);
	ProgressStep(1000, 1001, "%u significant anchors", g_HitCount);

	CloseStdioFile(g_fOut);
	CloseStdioFile(g_fTsv3Out);
	}
