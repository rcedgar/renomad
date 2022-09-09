#include "myutils.h"
#include "joinera.h"

static bool OnKmer(const string &KmerString,
  const vector<uint> &Counts, uint64 UserData)
	{
	JoinerA *J = (JoinerA *) UserData;
	J->OnKmer(KmerString, Counts);
	return true;
	}

void JoinerA::RunA(const vector<string> &KDBNames,
  Atab &A, fn_JoinerAOnAnchor OnAnchor, uint64 UserData)
	{
	m_OnAnchor = OnAnchor;
	m_UserData = UserData;
	m_A = &A;
	m_CurrentAnchor = "";
	Run(KDBNames, ::OnKmer, (uint64) this);
	if (m_A->m_TargetCount > 0)
		m_OnAnchor(*m_A, m_UserData);
	}

bool JoinerA::OnKmer(const string &KmerString,
  const vector<uint> &Counts)
	{
	asserta(SIZE(KmerString) == m_k);

	uint AnchorK = m_k/2;
	uint TargetK = m_k - AnchorK;

	string Anchor;
	for (uint i = 0; i < AnchorK; ++i)
		{
		char c = KmerString[i];
		Anchor += c;
		}

	string Target;
	for (uint i = 0; i < TargetK; ++i)
		{
		char c = KmerString[AnchorK + i];
		Target += c;
		}

	if (Anchor != m_CurrentAnchor)
		{
		if (m_A->m_TargetCount > 0)
			{
			bool Ok = m_OnAnchor(*m_A, m_UserData);
			if (!Ok)
				return false;
			}
		m_CurrentAnchor = Anchor;
		m_A->SetAnchor(Anchor.c_str());
		}

	m_A->AddTarget(Target.c_str(), Counts.data());

	return true;
	}
