#include "myutils.h"
#include "params.h"
#include <time.h>

static void Make_cs(uint SampleCount, vector<int> &cs)
	{
	cs.clear();
	for (uint i = 0; i < SampleCount; ++i)
		cs.push_back(i < SampleCount/2 ? 1 : -1);
	random_shuffle(cs.begin(), cs.end());
	}

static void Make_fs(uint TargetCount, vector<bool> &fs)
	{
	fs.clear();
	for (uint i = 0; i < TargetCount; ++i)
		fs.push_back(i < TargetCount/2);
	random_shuffle(fs.begin(), fs.end());
	}

double CalcPvalueBrute1(const vector<vector<uint> > &CountsVec,
  const vector<int> &c_js, const vector<bool> &fs, bool Trace)
	{
	const uint TargetCount = SIZE(CountsVec);
	asserta(TargetCount > 0);
	const uint SampleCount = SIZE(CountsVec[0]);
	asserta(SIZE(c_js) == SampleCount);
	asserta(SIZE(fs) == TargetCount);

	vector<uint> n_js;
	n_js.reserve(SampleCount);
	for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		uint n_j = 0;
		for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
			{
			uint n = CountsVec[TargetIndex][SampleIndex];
			n_j += n;
			}
		n_js.push_back(n_j);
		}

	if (Trace)
		{
		Log("n_js: ");
		for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			uint n_j = n_js[SampleIndex];
			Log(" %u", n_j);
			}
		Log("\n");
		}

	vector<double> mu_js;
	mu_js.reserve(SampleCount);
	uint M = 0;
	uint m = 0;
	for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		uint M_j = 0;
		uint m_j = 0;
		for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
			{
			uint n = CountsVec[TargetIndex][SampleIndex];
			M_j += n;
			if (fs[TargetIndex])
				m_j += n;
			}
		double mu_j = 0;
		if (M_j > 0)
			mu_j = double(m_j)/double(M_j);
		mu_js.push_back(mu_j);

		asserta(M_j == n_js[SampleIndex]);
		M += M_j;
		m += m_j;
		}
	double mu = 0;
	if (M > 0)
		mu = double(m)/double(M);

	if (Trace)
		{
		Log("mu_js: ");
		for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			double mu_j = mu_js[SampleIndex];
			Log(" %.2g", mu_j);
			}
		Log("\n");
		Log("mu = %.2g\n", mu);
		}

	vector<double> S_js;
	S_js.reserve(SampleCount);
	for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		uint n_j = n_js[SampleIndex];
		double mu_j = mu_js[SampleIndex];
		double S_j = sqrt(n_j)*(mu_j - mu);
		S_js.push_back(S_j);
		}

	if (Trace)
		{
		Log("S_js: ");
		for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			double S_j = S_js[SampleIndex];
			Log(" %.2g", S_j);
			}
		Log("\n");
		}

	double S = 0;
	double Sum_cj2 = 0;
	double Sum_cj2_ngt0 = 0;
	double Sum_cj_sqrtnj = 0;

	for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		uint n_j = n_js[SampleIndex];
		int c_j = c_js[SampleIndex];

		S += S_js[SampleIndex]*c_js[SampleIndex];
		Sum_cj2 += c_j*c_j;
		Sum_cj_sqrtnj += c_j*sqrt(n_j);
		if (n_j > 0)
			Sum_cj2_ngt0 += c_j*c_j;
		}
	asserta(Sum_cj2 == SampleCount);
	asserta(Sum_cj2_ngt0 == SampleCount);
	double P = 1;
	if (Sum_cj_sqrtnj == 0)
		{
		// double a = 0;
		double exp1_top = 2.0*S*S;
		double exp1 = exp1_top/Sum_cj2_ngt0;
		P = 2*exp(-exp1);
		}
	else
		{
		double ablob2 = float(M)*Sum_cj2/float(Sum_cj_sqrtnj*Sum_cj_sqrtnj);
		double ablob = sqrt(ablob2);
		double a = 1.0/(1 + ablob);
		double exp1_top = 2.0*(1 - a*a)*S*S;
		double exp1 = exp1_top/Sum_cj2_ngt0;
		double exp2_top = 2*a*a*M*S*S;
		double exp2 = exp2_top/(Sum_cj_sqrtnj*Sum_cj_sqrtnj);

		P = 2*exp(-exp1) + 2*exp(-exp2);
		}
	asserta(!isnan(P));
	if (P > 1)
		P = 1;
	if (Trace)
		Log("P = %.2g\n", P);
	return P;
	}

double CalcPvalueBrute(const vector<vector<uint> > &CountsVec)
	{
	const bool Trace = opt(trace_brute);
	uint SplitCount = g_SplitCount;
	uint FuncCount = g_FuncCount;

	const uint TargetCount = SIZE(CountsVec);
	if (TargetCount <= 1)
		return 1.0;

	const uint SampleCount = SIZE(CountsVec[0]);
	double MinP = 1.0;
	for (uint iL = 0; iL < SplitCount; ++iL)
		{
		vector<int> cs;
		Make_cs(SampleCount, cs);
		for (uint iK = 0; iK < FuncCount; ++iK)
			{
			vector<bool> fs;
			Make_fs(TargetCount, fs);

			double P = CalcPvalueBrute1(CountsVec, cs, fs, Trace);
			if (P < MinP)
				MinP = P;
			}
		}
	uint TestCount = SplitCount*FuncCount;
	double BonferroniCorrectedP = MinP*TestCount;
	if (BonferroniCorrectedP > 1)
		BonferroniCorrectedP = 1;
	return BonferroniCorrectedP;
	}

#if 1
static vector<vector<uint> > g_Table;

static void AddRow(uint j, uint n1, uint n2, uint n3, uint n4)
	{
	g_Table[0][j] = n1;
	g_Table[1][j] = n2;
	g_Table[2][j] = n3;
	g_Table[3][j] = n4;
	}

static void InitTable()
	{
	g_Table.resize(4);
	for (uint i = 0; i < 4; ++i)
		g_Table[i].resize(5);

	AddRow(0, 45, 5, 2, 1);
	AddRow(1, 93, 11, 3, 3);
	AddRow(2, 50, 6, 1, 2);
	AddRow(3, 2, 38, 1, 1);
	AddRow(4, 20, 370, 7, 8);
	}

static void TestBrute1()
	{
	InitTable();

	vector<bool> fs;
	fs.push_back(1);
	fs.push_back(0);
	fs.push_back(1);
	fs.push_back(0);

	vector<int> cs;
	cs.push_back(1);
	cs.push_back(1);
	cs.push_back(1);
	cs.push_back(-1);
	cs.push_back(-1);

	CalcPvalueBrute1(g_Table, cs, fs, true);
	}

// L is number of splits (cs), K is number of functions (fs)
static void TestBrute(uint L, uint K)
	{
	InitTable();

	const uint TargetCount = SIZE(g_Table);
	asserta(TargetCount > 0);
	const uint SampleCount = SIZE(g_Table[0]);

	srand(unsigned(time(0)));
	
	double MinP = 1.0;
	for (uint iL = 0; iL < L; ++iL)
		{
		vector<int> cs;
		Make_cs(SampleCount, cs);
		for (uint iK = 0; iK < K; ++iK)
			{
			vector<bool> fs;
			Make_fs(TargetCount, fs);

			double P = CalcPvalueBrute1(g_Table, cs, fs, false);
			if (P < MinP)
				MinP = P;
			}
		}
	uint TestCount = L*K;
	Log("\n");
	Log("MinP %.3g, corrected %.3g\n", MinP, MinP*TestCount);
	}

void cmd_test()
	{
	const string &NotUsed = opt(test);
//	TestBrute1();
	TestBrute(50, 10);

	//double REP = CalcREPvalue(g_Table);
	//ProgressLog("REP = %.3g\n", REP);
	}

#endif
