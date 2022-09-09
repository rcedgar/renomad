#include "myutils.h"
#include "joinercv.h"

static bool OnKmer(const string &KmerString,
  const vector<uint> &Counts, uint64 UserData)
	{
	JoinerCV *J = (JoinerCV *) UserData;
	J->OnKmer(KmerString, Counts);
	return true;
	}

void JoinerCV::RunCV(const vector<string> &KDBNames,
  fn_JoinerCVOnAnchor OnAnchor, uint64 UserData)
	{
	m_OnAnchor = OnAnchor;
	m_UserData = UserData;
	m_CurrentAnchor = "";
	Run(KDBNames, ::OnKmer, (uint64) this);
	if (!m_CountsVec.empty())
		m_OnAnchor(m_CurrentAnchor, m_Targets, m_CountsVec, m_UserData);
	}

bool JoinerCV::OnKmer(const string &KmerString,
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
		if (!m_CountsVec.empty())
			{
			bool Ok = m_OnAnchor(m_CurrentAnchor, m_Targets,
			  m_CountsVec, m_UserData);
			if (!Ok)
				return false;
			}
		m_CurrentAnchor = Anchor;
		m_Targets.clear();
		m_CountsVec.clear();
		}

	m_Targets.push_back(Target);
	m_CountsVec.push_back(Counts);

	return true;
	}
