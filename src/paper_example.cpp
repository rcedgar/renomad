#include "myutils.h"
#include "nomad.h"

double CalcPvalueBrute1(const vector<vector<uint> > &CountsVec,
  const vector<int> &c_js, const vector<bool> &fs, bool Trace);
void GetRandomKmer(uint k, string &Kmer);

static void AddRow(vector<vector<uint> > &CountsVec,
  uint j, uint n1, uint n2, uint n3, uint n4)
	{
	CountsVec[0][j] = n1;
	CountsVec[1][j] = n2;
	CountsVec[2][j] = n3;
	CountsVec[3][j] = n4;
	}

static void SetCountsVec(vector<vector<uint> > &CountsVec)
	{
	CountsVec.clear();
	CountsVec.resize(4);
	for (uint i = 0; i < 4; ++i)
		CountsVec[i].resize(5);

	AddRow(CountsVec, 0, 45, 5, 2, 1);
	AddRow(CountsVec, 1, 93, 11, 3, 3);
	AddRow(CountsVec, 2, 50, 6, 1, 2);
	AddRow(CountsVec, 3, 2, 38, 1, 1);
	AddRow(CountsVec, 4, 20, 370, 7, 8);
	}

void GetPaperExampleData(vector<vector<uint> > &CountsVec,
  vector<bool> &fs, vector<int> &cs)
	{
	SetCountsVec(CountsVec);

	fs.clear();
	fs.push_back(1);
	fs.push_back(0);
	fs.push_back(1);
	fs.push_back(0);

	cs.clear();
	cs.push_back(1);
	cs.push_back(1);
	cs.push_back(1);
	cs.push_back(-1);
	cs.push_back(-1);
	}

static void PaperExample_Brute()
	{
	vector<vector<uint> > CountsVec;
	vector<bool> fs;
	vector<int> cs;
	GetPaperExampleData(CountsVec, fs, cs);
	double P = CalcPvalueBrute1(CountsVec, cs, fs, true);

	Log("P_brute = %.3g\n", P);
	}


static void PaperExample_Nomad()
	{
	vector<vector<uint> > CountsVec;
	vector<bool> fs;
	vector<int> cs;
	GetPaperExampleData(CountsVec, fs, cs);

	vector<string> SampleNames;
	SampleNames.push_back("Sam1");
	SampleNames.push_back("Sam2");
	SampleNames.push_back("Sam3");
	SampleNames.push_back("Sam4");
	SampleNames.push_back("Sam5");

	uint TargetCount = 4;
	uint SampleCount = 5;

	Nomad Nd;
	Nd.Init(27, 27, SampleNames, 1000);

	string Anchor;
	vector<string> Targets;
	GetRandomKmer(27, Anchor);
	for (uint i = 0; i < SIZE(fs); ++i)
		{
		string Target;
		GetRandomKmer(27, Target);
		Targets.push_back(Target);
		}

	Nd.FromVecs(Anchor, Targets, CountsVec);

	Nd.m_SelectedSampleCount = SampleCount;
	Nd.m_SelectedTargetCount = TargetCount;

	for (uint i = 0; i < SampleCount; ++i)
		Nd.m_SelectedSamples[i] = i;

	for (uint i = 0; i < TargetCount; ++i)
		Nd.m_SelectedTargets[i] = i;

	Nd.m_fs = Nd.m_f_vec;
	for (uint i = 0; i < TargetCount; ++i)
		Nd.m_fs[i] = fs[i];

	Nd.m_cs = Nd.m_c_vec;
	for (uint i = 0; i < SampleCount; ++i)
		Nd.m_cs[i] = cs[i];

	Nd.Set_njs();
	Nd.Set_mujs();
	Nd.Set_Sjs();
	double P1 = Nd.CalcP1();

	Log("P_nomad1 = %.3g\n", P1);

	double P = Nd.CalcP();
	Log("P_nomad = %.3g\n", P);
	}

#if 0
void cmd_test()
	{
	opt(test);
	PaperExample_Brute();
	PaperExample_Nomad();
	}
#endif // 0
