#include "myutils.h"
#include "joiner.h"

void Joiner::Clear()
	{
	const uint N = SIZE(m_KDBs);
	for (uint i = 0; i < N; ++i)
		{
		CKMCFile *KDB = m_KDBs[i];
		asserta(KDB != 0);
		KDB->Close();
		delete KDB;
		}
	m_KDBs.clear();
	m_KDBNames.clear();
	m_MinCount = 1;
	m_MaxCount = UINT_MAX;
	m_k = UINT_MAX;
	if (m_Kmer != 0)
		{
		delete m_Kmer;
		m_Kmer = 0;
		}
	}

void Joiner::Open(const vector<string> &KDBNames)
	{
	Clear();

	m_KDBNames = KDBNames;
	const uint N = SIZE(m_KDBNames);
	asserta(N > 1);

	for (uint i = 0; i < N; ++i)
		{
		CKMCFile *KDB = new CKMCFile;
		const string &KDBName = m_KDBNames[i];

		if (m_MinCount != 1)
			{
			bool Ok = KDB->SetMinCount(m_MinCount);
			asserta(Ok);
			}

		if (m_MaxCount != UINT_MAX)
			{
			bool Ok = KDB->SetMaxCount(m_MaxCount);
			asserta(Ok);
			}

		bool Ok = KDB->OpenForListing(KDBName);
		if (!Ok)
			Die("KDB->OpenForListing(%s) = false", KDBName.c_str());

		uint32 _kmer_length;
		uint32 _mode;
		uint32 _counter_size;
		uint32 _lut_prefix_length;
		uint32 _signature_len;
		uint32 _min_count;
		uint64 _max_count;
		uint64 _total_kmers;

		KDB->Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, 
		 _signature_len, _min_count, _max_count, _total_kmers);
		if (m_KDBs.empty())
			m_k = _kmer_length;
		else
			asserta(_kmer_length == m_k);

		m_KDBs.push_back(KDB);

		m_CurrentKmers.push_back("");
		m_CurrentCounts.push_back(0);
		m_TotalKmers.push_back(_total_kmers);
		m_ScannedKmers.push_back(0);
		}
	m_Kmer = new CKmerAPI(m_k);

	m_CurrentKmer.clear();
	for (uint i = 0; i < N; ++i)
		{
		GetNext(i);
		const string &ks = m_CurrentKmers[i];
		if (m_CurrentKmer.empty() || ks < m_CurrentKmer)
			m_CurrentKmer = ks;
		}
	Validate();
	}

bool Joiner::GetNext(uint KDBIndex)
	{
	asserta(KDBIndex < SIZE(m_KDBs));
	CKMCFile *KDB = m_KDBs[KDBIndex];
	uint64 Count64;
	bool Ok = KDB->ReadNextKmer(*m_Kmer, Count64);
	if (!Ok)
		{
		m_CurrentKmers[KDBIndex] = "$";
		return false;
		}
	
	m_ScannedKmers[KDBIndex] += 1;
	m_CurrentCounts[KDBIndex] =
	  (Count64 > UINT_MAX ? UINT_MAX : uint32(Count64));

	string OldKmer = m_CurrentKmers[KDBIndex];
	m_Kmer->to_string(m_CurrentKmers[KDBIndex]);
	const string &NewKmer = m_CurrentKmers[KDBIndex];

	if (!OldKmer.empty() && !NewKmer.empty())
		{
		if (NewKmer <= OldKmer)
			{
			Log("\n");
			Log("Name %s\n", m_KDBNames[KDBIndex].c_str());
			Log("%s  old\n", OldKmer.c_str());
			Log("%s  new\n", NewKmer.c_str());
			Die("GetNext(%u) New <= Old", KDBIndex);
			}
		}

	return true;
	}

void Joiner::GetCurrent(string &KmerString, vector<uint> &Counts) const
	{
	KmerString = m_CurrentKmer;

	Counts.clear();
	const uint N = SIZE(m_CurrentKmers);
	for (uint i = 0; i < N; ++i)
		{
		const string &ks = m_CurrentKmers[i];
		if (ks > m_CurrentKmer || ks == "$")
			Counts.push_back(0);
		else if (ks == m_CurrentKmer)
			Counts.push_back(m_CurrentCounts[i]);
		else
			asserta(false);
		}
	}

bool Joiner::Advance()
	{
	const uint N = SIZE(m_CurrentKmers);

	const string PrevCurrentKmer = m_CurrentKmer;
	for (uint i = 0; i < N; ++i)
		{
		const string &ks = m_CurrentKmers[i];
		if (ks == PrevCurrentKmer)
			GetNext(i);
		}

	string NextCurrentKmer;
	for (uint i = 0; i < N; ++i)
		{
		const string &ks = m_CurrentKmers[i];
		if (ks == "$")
			continue;
		if (NextCurrentKmer == "")
			NextCurrentKmer = ks;
		if (ks > PrevCurrentKmer && ks < NextCurrentKmer)
			NextCurrentKmer = ks;
		}
	if (NextCurrentKmer == "")
		return false;

	if (NextCurrentKmer != PrevCurrentKmer)
		asserta(NextCurrentKmer > PrevCurrentKmer);

	bool AnyMore = (NextCurrentKmer > PrevCurrentKmer);
	m_CurrentKmer = NextCurrentKmer;
	Validate();

	if (!AnyMore)
		return false;
	return true;
	}

void Joiner::Validate() const
	{
	const uint N = SIZE(m_CurrentKmers);
	asserta(SIZE(m_KDBs) == N);
	asserta(SIZE(m_CurrentKmers) == N);
	asserta(SIZE(m_CurrentCounts) == N);
	asserta(!m_CurrentKmer.empty());

	for (uint i = 0; i < N; ++i)
		{
		const string &ks = m_CurrentKmers[i];
		if (ks != "$" && ks < m_CurrentKmer)
			{
			Log("\n");
			Log("%s  [%u]\n", ks.c_str(), i);
			Log("%s  current\n", m_CurrentKmer.c_str());
			Die("ks < m_CurrentKmer");
			}
		}
	}

double Joiner::GetPctDone() const
	{
	const uint N = SIZE(m_CurrentKmers);
	double SumPct = 0;
	for (uint i = 0; i < N; ++i)
		{
		uint64 TotalKmers = m_TotalKmers[i];
		uint64 ScannedKmers = m_ScannedKmers[i];
		asserta(ScannedKmers <= TotalKmers);
		double Pct = GetPct(double(ScannedKmers), double(TotalKmers));
		SumPct += Pct;
		}
	double AvgPct = SumPct/N;
	return AvgPct;
	}

void Joiner::Run(const vector<string> &KDBNames,
  fn_JoinerOnKmer OnKmer, uint64 UserData)
	{
	m_OnKmer = OnKmer;
	m_UserData = UserData;

	const uint N = SIZE(KDBNames);

	Open(KDBNames);

	string KmerString;
	vector<uint> Counts;
	for (;;)
		{
		GetCurrent(KmerString, Counts);
		bool Ok1 = Advance();
		if (!Ok1)
			break;
		bool Ok2 = m_OnKmer(KmerString, Counts, m_UserData);
		if (!Ok2)
			break;
		}

	Clear();
	}
