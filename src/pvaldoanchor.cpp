#include "myutils.h"
#include "params.h"
#include "nomad.h"

double CalcPvalue(const vector<vector<uint> > &CountsVec);
void MxOut(FILE *f, const string &Anchor, 
  const vector<string> &TargetVec, 
  const vector<vector<uint> > &CountsVec);

static Nomad g_Nomad;
static bool g_InitDone = false;

void InitNomad(uint ka, uint kt, const vector<string> &SampleNames)
	{
	g_Nomad.Init(ka, kt, SampleNames, 1000);
	}

void PvalDoAnchor(FILE *fOut, FILE *fMxOut, 
  const string &Anchor, 
  const vector<string> &TargetVec, 
  const vector<vector<uint> > &CountsVec)
	{
//	double P = CalcPvalue(CountsVec);

	g_Nomad.FromVecs(Anchor, TargetVec, CountsVec);
	double P2 = g_Nomad.CalcP();

//	Log("%.3g %.3g\n", P, P2);

	if (P2 <= g_MaxP)
		{
		Pf(fOut, "%.2g\t%s\n", P2, Anchor.c_str());
		MxOut(fMxOut, Anchor, TargetVec, CountsVec);
		}
	}
