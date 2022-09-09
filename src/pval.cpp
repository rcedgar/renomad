#include "myutils.h"
#include "params.h"

static FILE *g_fOut;
static FILE *g_fMxOut;

void InitNomad(uint ka, uint kt, const vector<string> &SampleNames);

double CalcPvalueBrute(const vector<vector<uint> > &CountsVec);

void PvalDoAnchor(FILE *fOut, FILE *fMx, 
  const string &Anchor, 
  const vector<string> &TargetVec, 
  const vector<vector<uint> > &CountsVec);

static uint GetAnchorTotal(const vector<vector<uint> > &CountsVec)
	{
	const uint TargetCount = SIZE(CountsVec);
	if (TargetCount == 0)
		return 0;

	const uint SampleCount = SIZE(CountsVec[0]);
	uint Sum = 0;
	for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			uint Count = CountsVec[TargetIndex][SampleIndex];
			Sum += Count;
			}
		}
	return Sum;
	}

static void RemoveLowCountTargets(const vector<vector<uint> > &CountsVecIn,
  vector<vector<uint> > &CountsVecOut)
	{
	CountsVecOut.clear();
	const uint TargetCount = SIZE(CountsVecIn);
	if (TargetCount <= 1)
		return;

	const uint SampleCount = SIZE(CountsVecIn[0]);
	for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		uint Total = 0;
		const vector<uint> &Counts = CountsVecIn[TargetIndex];
		asserta(SIZE(Counts) == SampleCount);
		for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			uint Count = CountsVecIn[TargetIndex][SampleIndex];
			Total += Count;
			}
		if (Total >= g_MinTargetTotal)
			CountsVecOut.push_back(Counts);
		}
	}

static void RemoveLowCountSamples(const vector<vector<uint> > &CountsVecIn,
  vector<vector<uint> > &CountsVecOut)
	{
	CountsVecOut.clear();
	const uint TargetCount = SIZE(CountsVecIn);
	if (TargetCount <= 1)
		return;

	const uint SampleCount = SIZE(CountsVecIn[0]);
	vector<bool> Discards;
	uint KeepCount = 0;
	for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		uint SampleTotal = 0;
		for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
			{
			uint Count = CountsVecIn[TargetIndex][SampleIndex];
			SampleTotal += Count;
			}
		bool Discard = (SampleTotal < g_MinSampleTotal);
		Discards.push_back(Discard);
		if (!Discard)
			++KeepCount;
		}

	if (KeepCount < g_MinSamplesPerAnchor)
		return;

	CountsVecOut.resize(TargetCount);
	for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		const vector<uint> &CountsIn = CountsVecIn[TargetIndex];
		vector<uint> &CountsOut = CountsVecOut[TargetIndex];
		for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			if (Discards[SampleIndex])
				continue;
			uint Count = CountsIn[SampleIndex];
			CountsOut.push_back(Count);
			}
		}
	}

void MxOut(FILE *f, const string &Anchor, 
  const vector<string> &TargetVec, 
  const vector<vector<uint> > &CountsVec)
	{
	if (f == 0)
		return;

	const uint TargetCount = SIZE(TargetVec);
	asserta(SIZE(CountsVec) == TargetCount);
	for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		fprintf(f, "%s", Anchor.c_str());
		fprintf(f, "\t%s", TargetVec[TargetIndex].c_str());

		const vector<uint> &Counts = CountsVec[TargetIndex];
		const uint SampleCount = SIZE(Counts);
		for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			uint Count = Counts[SampleIndex];
			fprintf(f, "\t%u", Count);
			}
		fprintf(f, "\n");
		}
	}

static void Prune(const vector<vector<uint> > &CountsVecIn,
  vector<vector<uint> > &CountsVecOut)
	{
	vector<vector<uint> > CountsVecPruned1;
	RemoveLowCountSamples(CountsVecIn, CountsVecPruned1);
	RemoveLowCountTargets(CountsVecPruned1, CountsVecOut);
	}

double CalcPvalue(const vector<vector<uint> > &CountsVec)
	{
	uint TargetCountBeforeFiltering = SIZE(CountsVec);
	if (TargetCountBeforeFiltering < g_MinTargetsPerAnchor)
		return 1.0;

	uint AnchorTotalBeforeFiltering = GetAnchorTotal(CountsVec);
	if (AnchorTotalBeforeFiltering < g_MinAnchorTotalBeforeFiltering)
		return 1.0;

	vector<vector<uint> > CountsVecPruned;
	Prune(CountsVec, CountsVecPruned);
	uint TargetCountAfterFiltering = SIZE(CountsVecPruned);
	if (TargetCountAfterFiltering < g_MinTargetsPerAnchor)
		return 1.0;

	double P = CalcPvalueBrute(CountsVecPruned);
	return P;
	}

void cmd_pval()
	{
	const string &InputFileName = opt(pval);
	const string &OutputFileName = opt(output);
	const string &MxOutputFileName = opt(mxout);

	SetParams_CmdLine();

	g_fOut = CreateStdioFile(OutputFileName);
	g_fMxOut = CreateStdioFile(MxOutputFileName);

	FILE *fIn = OpenStdioFile(InputFileName);
	string Line;
	bool Ok = ReadLineStdioFile(fIn, Line);
	asserta(Ok);

	vector<string> HdrFields;
	Split(Line, HdrFields, '\t');
	const uint HdrFieldCount = SIZE(HdrFields);
	asserta(HdrFieldCount > 2);
	asserta(HdrFields[0] == "Anchor");
	asserta(HdrFields[1] == "Target");
	const uint SampleCount = HdrFieldCount - 2;

	vector<string> SampleNames;
	for (uint j = 0; j < SampleCount; ++j)
		SampleNames.push_back(HdrFields[j+2]);
	InitNomad(27, 27, SampleNames);

	vector<string> Fields;
	string CurrAnchor;
	vector<string> TargetVec;
	vector<vector<uint> > CountsVec;
	uint TargetCount = 0;
	while (ReadLineStdioFile(fIn, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == HdrFieldCount);

		const string &Anchor = Fields[0];
		const string &Target = Fields[1];

		if (Anchor != CurrAnchor)
			{
			if (!CurrAnchor.empty())
				PvalDoAnchor(g_fOut, g_fMxOut,
				  CurrAnchor, TargetVec, CountsVec);

			TargetVec.clear();
			CountsVec.clear();
			CurrAnchor = Anchor;
			TargetCount = 0;
			}

		vector<uint> Counts;
		uint TargetMaxCount = 0;
		uint TargetTotal = 0;
		for (uint i = 0; i < SampleCount; ++i)
			{
			uint Count = StrToUint(Fields[i+2]);
			if (Count > TargetMaxCount)
				TargetMaxCount = Count;
			Counts.push_back(Count);
			}

		if (TargetMaxCount >= g_MinTargetMaxCount)
			{
			TargetVec.push_back(Target);
			CountsVec.push_back(Counts);
			TargetCount += 1;
			}
		}
	PvalDoAnchor(g_fOut, g_fMxOut,
	  CurrAnchor, TargetVec, CountsVec);

	CloseStdioFile(fIn);
	CloseStdioFile(g_fOut);
	CloseStdioFile(g_fMxOut);
	}
