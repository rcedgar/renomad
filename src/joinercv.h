#pragma once

#include "joiner.h"

typedef bool fn_JoinerCVOnAnchor(
  const string &AnchorString,
  const vector<string> &Targets,
  const vector<vector<uint> > &CountsVec,
  uint64 UserData);

class JoinerCV : public Joiner
	{
public:
	string m_CurrentAnchor;
	vector<string> m_Targets;
	vector<vector<uint> > m_CountsVec;
	fn_JoinerCVOnAnchor *m_OnAnchor = 0;
	uint64 m_UserData;

public:
	void RunCV(const vector<string> &KDBNames,
	  fn_JoinerCVOnAnchor OnAnchor, uint64 UserData);

	bool OnKmer(const string &KmerString,
	  const vector<uint> &Counts);
	};
