#include "myutils.h"
#include "joiner.h"

void ReadStringsFromFile(const string &FileName,
  vector<string> &Strs, bool DeleteEmptyStrings)
	{
	Strs.clear();
	FILE *f = OpenStdioFile(FileName);
	string s;
	while (ReadLineStdioFile(f, s))
		Strs.push_back(s);
	CloseStdioFile(f);
	}

void cmd_join()
	{
	const string &InputFileName = opt(join);
	const string &OutputFileName = opt(output);

	uint MinSum = 32;
	uint MinMaxCount = 4;

	vector<string> KDBNames;
	ReadStringsFromFile(InputFileName, KDBNames, true);
	const uint N = SIZE(KDBNames);

	Joiner J;
	J.Open(KDBNames);

	string KmerString;
	vector<uint> Counts;
	uint64 ProcessedCount = 0;
	uint64 OutputCount = 0;
	FILE *fOut = CreateStdioFile(OutputFileName);
	uint64 BytesWritten = 0;
	for (;;)
		{
		if (ProcessedCount%1000000 == 0)
			{
			double Pct = GetPct(double(OutputCount), double(ProcessedCount));
			const char *ns = Int64ToStr(ProcessedCount);
			double PctDone = J.GetPctDone();
			string BytesWrittenStr = string(MemBytesToStr(BytesWritten));

			uint64 EstTotal = PctDone == 0 ? 0 : uint64(BytesWritten*100.0/PctDone);
			string EstTotalStr = string(MemBytesToStr(EstTotal));
			Progress("%.2f%% done, processed %s, output %.1f%% (%s / %s)\r",
			  PctDone, ns, Pct, BytesWrittenStr.c_str(), EstTotalStr.c_str());
			}

		J.GetCurrent(KmerString, Counts);
		uint Sum = 0;
		uint MaxCount = 0;
		for (uint i = 0; i < N; ++i)
			{
			uint n = Counts[i];
			if (n > MaxCount)
				MaxCount = n;
			Sum += n;
			}
		if (Sum >= MinSum && MaxCount >= MinMaxCount)
			{
			if (fOut != 0)
				{
				++OutputCount;

				string OutStr = KmerString;
				for (uint i = 0; i < N; ++i)
					Psa(OutStr, "\t%u", Counts[i]);
				OutStr += '\n';
				fputs(OutStr.c_str(), fOut);
				BytesWritten += SIZE(OutStr);
				}
			}
		bool Ok = J.Advance();
		if (!Ok)
			break;
		++ProcessedCount;
		}

	CloseStdioFile(fOut);
	J.Clear();
	}
