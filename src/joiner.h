#pragma once

#include "kmc_file.h"

typedef bool fn_JoinerOnKmer(const string &KmerString,
  const vector<uint> &Counts, uint64 UserData);

class Joiner
	{
public:
	vector<CKMCFile *> m_KDBs;
	vector<string> m_KDBNames;
	uint m_MinCount = 1;
	uint m_MaxCount = UINT_MAX;
	uint m_k = UINT_MAX;

	CKmerAPI *m_Kmer = 0;

	string m_CurrentKmer;
	vector<string> m_CurrentKmers;
	vector<uint> m_CurrentCounts;
	vector<uint64> m_TotalKmers;
	vector<uint64> m_ScannedKmers;

	fn_JoinerOnKmer *m_OnKmer = 0;
	uint64 m_UserData = 0;

public:
	void Run(const vector<string> &KDBNames,
	  fn_JoinerOnKmer OnKmer, uint64 UserData);

	void Clear();
	void Open(const vector<string> &KDBNames);
	void GetCurrent(string &KmerString,
	  vector<uint> &Counts) const;
	bool Advance();
	void Validate() const;
	double GetPctDone() const;

private:
	bool GetNext(uint KDBIndex);
	};
