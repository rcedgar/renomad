#pragma once

#include "joiner.h"
#include "atab.h"

typedef bool fn_JoinerAOnAnchor(Atab &A, uint64 UserData);

class JoinerA : public Joiner
	{
public:
	Atab *m_A = 0;
	fn_JoinerAOnAnchor *m_OnAnchor = 0;
	uint64 m_UserData = 0;
	string m_CurrentAnchor;

public:
	void RunA(const vector<string> &KDBNames,
	  Atab &A, fn_JoinerAOnAnchor OnAnchor,
	  uint64 UserData);

	bool OnKmer(const string &KmerString,
	  const vector<uint> &Counts);
	};
