#include "myutils.h"
#include "atab.h"
#include "params.h"

void GetRandomKmer(uint k, string &Kmer)
	{
	Kmer.clear();
	for (uint i = 0; i < k; ++i)
		{
		uint r = randu32();
		char c = "ACGT"[r%4];
		Kmer.push_back(c);
		}
	}

void Atab::Init(uint ka, uint kt, 
 const vector<string> &SampleNames,
 uint InitialMaxTargetCount)
	{
	Clear();

	m_SampleNames = SampleNames;
	m_SampleCount = SIZE(m_SampleNames);
	m_ka = ka;
	m_kt = kt;
	m_TargetCount = 0;
	m_MaxTargetCount = InitialMaxTargetCount;

	m_Anchor = myalloc(char, m_ka+1);
	m_Anchor[m_ka] = 0;

	m_Targets = myalloc(char *, m_MaxTargetCount);
	m_Counts = myalloc(uint *, m_MaxTargetCount);
	for (uint i = 0; i < m_MaxTargetCount; ++i)
		{
		m_Targets[i] = myalloc(char, m_kt+1);
		m_Targets[i][kt] = 0;

		m_Counts[i] = myalloc(uint, m_SampleCount);
		}

	m_TargetTotals = myalloc(uint, m_MaxTargetCount);
	m_SampleTotals = myalloc(uint, m_SampleCount);
	}

void Atab::SetAnchor(const char *Anchor)
	{
	memcpy(m_Anchor, Anchor, m_ka);
	m_TargetCount = 0;
	m_Total = 0;
	zero(m_SampleTotals, m_SampleCount);
	}

void Atab::AddTarget(const char *Target, const uint *Counts)
	{
	uint Sum = 0;
	uint MaxCount = 0;
	for (uint j = 0; j < m_SampleCount; ++j)
		{
		uint n = Counts[j];
		if (n > MaxCount)
			MaxCount = n;
		Sum += n;
		}
	if (Sum == 0)
		return;
	if (Sum < g_MinTargetTotal)
		return;
	if (MaxCount < g_MinTargetMaxCount)
		return;

	if (m_TargetCount == m_MaxTargetCount)
		{
		Expand();
		asserta(m_TargetCount < m_MaxTargetCount);
		}

	m_Total += Sum;
	m_TargetTotals[m_TargetCount] = Sum;
	for (uint j = 0; j < m_SampleCount; ++j)
		{
		uint n = Counts[j];
		m_SampleTotals[j] += n;
		}

	size_t Bytes = m_SampleCount*sizeof(Counts[0]);
	memcpy(m_Counts[m_TargetCount], Counts, Bytes);
	memcpy(m_Targets[m_TargetCount], Target, m_kt);
	++m_TargetCount;
	}

void Atab::Expand()
	{
	uint NewMaxTargetCount = m_MaxTargetCount + 256;

	uint *NewTargetTotals = myalloc(uint, NewMaxTargetCount);

	memcpy(NewTargetTotals, m_TargetTotals, m_TargetCount*sizeof(uint));
	myfree(m_TargetTotals);
	m_TargetTotals = NewTargetTotals;

	char **NewTargets = myalloc(char *, NewMaxTargetCount);
	uint **NewCounts = myalloc(uint *, NewMaxTargetCount);
	for (uint i = 0; i < NewMaxTargetCount; ++i)
		{
		NewCounts[i] = myalloc(uint, m_SampleCount);
		NewTargets[i] = myalloc(char, m_kt + 1);
		NewTargets[i][m_kt] = 0;
		}

	for (uint i = 0; i < m_TargetCount; ++i)
		{
		memcpy(NewTargets[i], m_Targets[i], m_kt + 1);
		memcpy(NewCounts[i], m_Counts[i], m_SampleCount*sizeof(uint));
		}

	for (uint i = 0; i < m_MaxTargetCount; ++i)
		{
		myfree(m_Targets[i]);
		myfree(m_Counts[i]);
		}
	myfree(m_Targets);
	myfree(m_Counts);

	m_Targets = NewTargets;
	m_Counts = NewCounts;

	m_MaxTargetCount = NewMaxTargetCount;
	}

void Atab::ToFileMx(FILE *f) const
	{
	if (f == 0)
		return;

	for (uint TargetIndex = 0; TargetIndex < m_TargetCount; ++TargetIndex)
		{
		fputs(m_Anchor, f);
		fputc('\t', f);
		fputs(m_Targets[TargetIndex], f);
		const uint *ns = m_Counts[TargetIndex];
		for (uint SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
			{
			uint n = ns[SampleIndex];
			fprintf(f, "\t%u", n);
			}
		fputc('\n', f);
		}
	}

void Atab::ToFileTsv3(FILE *f) const
	{
	if (f == 0)
		return;

	for (uint TargetIndex = 0; TargetIndex < m_TargetCount; ++TargetIndex)
		{
		const uint *ns = m_Counts[TargetIndex];
		for (uint SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
			{
			uint n = ns[SampleIndex];
			if (n > 0)
				{
				const char *SampleName = m_SampleNames[SampleIndex].c_str();

				fprintf(f, "%u", n);
				fputc('\t', f);
				fputs(m_Anchor, f);
				fputs(m_Targets[TargetIndex], f);
				fputc('\t', f);
				fputs(SampleName, f);
				fputc('\n', f);
				}
			}
		}
	}

void Atab::HdrToFileMx(FILE *f) const
	{
	if (f == 0)
		return;

	fputs("Anchor\tTarget", f);
	for (uint SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
		{
		const char *SampleName = m_SampleNames[SampleIndex].c_str();
		fputc('\t', f);
		fputs(SampleName, f);
		}
	fputc('\n', f);
	}

void Atab::FromVec(const vector<vector<uint> > &CountsVec)
	{
	string Anchor;
	GetRandomKmer(m_ka, Anchor);
	SetAnchor(Anchor.c_str());

	vector<string> Targets;
	uint TargetCount = SIZE(CountsVec);
	asserta(TargetCount <= m_MaxTargetCount);
	for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		string Target;
		GetRandomKmer(m_kt, Target);
		const vector<uint> &Counts = CountsVec[TargetIndex];
		AddTarget(Target.c_str(), Counts.data());
		}
	}

void Atab::FromVecs(const string &Anchor,
  const vector<string> &Targets,
  const vector<vector<uint> > &CountsVec)
	{
	asserta(SIZE(Anchor) == m_ka);
	SetAnchor(Anchor.c_str());

	uint TargetCount = SIZE(Targets);
	asserta(SIZE(CountsVec) == TargetCount);
	asserta(TargetCount <= m_MaxTargetCount);
	for (uint TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		const string &Target = Targets[TargetIndex];
		asserta(SIZE(Target) == m_kt);
		const vector<uint> &Counts = CountsVec[TargetIndex];
		AddTarget(Target.c_str(), Counts.data());
		}
	}
