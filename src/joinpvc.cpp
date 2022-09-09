#include "myutils.h"
#include "joinercv.h"

void ReadStringsFromFile(const string &FileName,
  vector<string> &Strs, bool DeleteEmptyStrings);
double CalcPvalueBrute(const vector<vector<uint> > &CountsVec);

static uint g_HitCount = 0;
static FILE *g_fOut = 0;
static double g_LastPctDone = 0;

static bool OnAnchor(
  const string &Anchor,
  const vector<string> &Targets,
  const vector<vector<uint> > &CountsVec,
  uint64 UserData)
	{
	double P = CalcPvalueBrute(CountsVec);
	if (P < 0.5)
		{
		JoinerCV *JC = (JoinerCV *) UserData;
		double Pct = JC->GetPctDone();
		if (Pct - g_LastPctDone > 0.1)
			{
			Progress("%u significant anchors, %.3g%% done\r",
			  ++g_HitCount, Pct);
			g_LastPctDone = Pct;
			}
//		Log("%8.2g  %s\n", P, Anchor.c_str());
		Pf(g_fOut, "%.2g\t%s\n", P, Anchor.c_str());
		}
	return true;
	}

void cmd_joinpvc()
	{
	const string &InputFileName = opt(joinpvc);
	const string &OutputFileName = opt(output);
	
	g_fOut = CreateStdioFile(OutputFileName);

	vector<string> KDBNames;
	ReadStringsFromFile(InputFileName, KDBNames, true);
	const uint N = SIZE(KDBNames);

	g_HitCount = 0;
	g_LastPctDone = 0;

	JoinerCV JC;
	JC.RunCV(KDBNames, OnAnchor, (uint64) &JC);

	CloseStdioFile(g_fOut);
	}
