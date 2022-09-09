#include "myutils.h"
#include "mxreader.h"
#include "params.h"

void MxReader::Run(const string &FileName, 
  Atab &A, uint ka, uint kt, 
  fn_MxReaderOnStart OnStart,
  fn_MxReaderOnAnchor OnAnchor,
  uint64 UserData)
	{
	m_UserData = UserData;
	m_fIn = OpenStdioFile(FileName);
	string Line;
	bool Ok = ReadLineStdioFile(m_fIn, Line);
	asserta(Ok);

	vector<string> HdrFields;
	Split(Line, HdrFields, '\t');
	const uint HdrFieldCount = SIZE(HdrFields);
	asserta(HdrFieldCount > 2);
	asserta(HdrFields[0] == "Anchor");
	asserta(HdrFields[1] == "Target");
	const uint SampleCount = HdrFieldCount - 2;
	asserta(SampleCount > 1);

	vector<string> SampleNames;
	for (uint i = 0; i < SampleCount; ++i)
		SampleNames.push_back(HdrFields[i+2]);

	m_A = &A;
	m_A->Init(ka, kt, SampleNames, g_InitialMaxTargetCount);
	if (OnStart != 0)
		{
		bool Ok = OnStart(ka, kt, SampleNames, UserData);
		if (!Ok)
			{
			CloseStdioFile(m_fIn);
			return;
			}
		}

	uint *Counts = myalloc(uint, SampleCount);

	vector<string> Fields;
	string CurrAnchor;
	uint TargetCount = 0;
	while (ReadLineStdioFile(m_fIn, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == HdrFieldCount);

		const string &Anchor = Fields[0];
		const string &Target = Fields[1];

		if (Anchor != CurrAnchor)
			{
			if (m_A->m_TargetCount > 0)
				{
				bool Ok = OnAnchor(*m_A, m_UserData);
				if (!Ok)
					{
					CloseStdioFile(m_fIn);
					return;
					}
				}

			m_A->SetAnchor(Anchor.c_str());
			CurrAnchor = Anchor;
			}

		uint TargetMaxCount = 0;
		for (uint i = 0; i < SampleCount; ++i)
			{
			uint Count = StrToUint(Fields[i+2]);
			if (Count > TargetMaxCount)
				TargetMaxCount = Count;
			Counts[i] = Count;
			}

		if (TargetMaxCount >= g_MinTargetMaxCount)
			m_A->AddTarget(Target.c_str(), Counts);
		}
	if (m_A->m_TargetCount > 0)
		OnAnchor(*m_A, m_UserData);

	CloseStdioFile(m_fIn);
	}

#if 0
static FILE *g_fOut;
static bool g_HdrDone = false;
static bool OnAnchor(const Atab &A, uint64 UserData)
	{
	if (!g_HdrDone)
		{
		A.HdrToFileMx(g_fOut);
		g_HdrDone = true;
		}
//	Log("A=%s\n", A.m_Anchor);
	A.ToFileMx(g_fOut);
	return true;
	}

void cmd_test()
	{
	string FileName = opt(test);
	g_fOut = CreateStdioFile(opt(output));

	MxReader MR;
	MR.Run(FileName, 27, 27, OnAnchor, 0);

	CloseStdioFile(g_fOut);
	}
#endif // 0
