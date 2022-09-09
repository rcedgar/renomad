#pragma once

#include "atab.h"

typedef bool fn_MxReaderOnStart(uint ka, uint kt,
  const vector<string> &SampleNames, uint64 UserData);

typedef bool fn_MxReaderOnAnchor(
  const Atab &A, uint64 UserData);

class MxReader
	{
public:
	string m_FileName;
	FILE *m_fIn = 0;
	uint64 m_UserData = 0;
	Atab *m_A = 0;

public:
	void Clear()
		{
		if (m_fIn != 0)
			CloseStdioFile(m_fIn);
		m_fIn = 0;
		m_UserData = 0;
		m_A = 0;
		}

	void Run(const string &FileName,
		Atab &A, uint ka, uint kt,
		fn_MxReaderOnStart OnStart,
		fn_MxReaderOnAnchor OnAnchor,
		uint64 UserData);
	};
