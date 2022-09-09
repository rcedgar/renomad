#pragma once

class Atab
	{
public:
	uint m_ka = 0;
	uint m_kt = 0;
	char *m_Anchor = 0;
	vector<string> m_SampleNames;
	uint m_SampleCount = 0;
	uint m_TargetCount = 0;
	uint m_MaxTargetCount = 0;
	char **m_Targets = 0;
	uint **m_Counts = 0;
	uint *m_SampleTotals = 0;
	uint *m_TargetTotals = 0;
	uint m_Total = 0;

public:
	void Clear()
		{
		myfree(m_Anchor);
		m_Anchor = 0;

		for (uint i = 0; i < m_MaxTargetCount; ++i)
			{
			myfree(m_Targets[i]);
			myfree(m_Counts[i]);
			}
		myfree(m_Targets);
		myfree(m_Counts);

		m_SampleNames.clear();
		m_Targets = 0;
		m_Counts = 0;
		m_MaxTargetCount = 0;

		m_ka = 0;
		m_kt = 0;
		}

	void Init(uint ka, uint kt, 
	  const vector<string> &SampleNames,
	  uint InitialMaxTargetCount);

	void SetAnchor(const char *Anchor);
	void AddTarget(const char *Target, const uint *Counts);
	void Expand();
	void HdrToFileMx(FILE *f) const;
	void ToFileMx(FILE *f) const;
	void ToFileTsv3(FILE *f) const;

	void FromVec(const vector<vector<uint> > &CountsVec);

	void FromVecs(const string &Anchor,
	  const vector<string> &Targets,
	  const vector<vector<uint> > &CountsVec);
	};
