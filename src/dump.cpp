#include "myutils.h"
#include "kmc_file.h"

void cmd_dump()
	{
	const string &input_file_name = opt(dump);
	const string &output_file_name = opt(output);

	CKMCFile kmer_data_base;
	uint32 min_count = 1;
	uint32 max_count = UINT_MAX;

	FILE *out_file = CreateStdioFile(output_file_name);

	bool Ok = kmer_data_base.OpenForListing(input_file_name);
	if (!Ok)
		Die("kmer_data_base.OpenForListing(%s)=false", input_file_name.c_str());

	uint32 _kmer_length;
	uint32 _mode;
	uint32 _counter_size;
	uint32 _lut_prefix_length;
	uint32 _signature_len;
	uint32 _min_count;
	uint64 _max_count;
	uint64 _total_kmers;

	kmer_data_base.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, 
	 _signature_len, _min_count, _max_count, _total_kmers);

	if (min_count != 1)
		{
		Ok = kmer_data_base.SetMinCount(min_count);
		asserta(Ok);
		}

	if (max_count != UINT_MAX)
		{
		Ok = kmer_data_base.SetMaxCount(max_count);
		asserta(Ok);
		}

	uint64 counter;
	char str[1024];
	CKmerAPI kmer_object(_kmer_length);
	string PrevKmer;
	while (kmer_data_base.ReadNextKmer(kmer_object, counter))
		{
		kmer_object.to_string(str);
		str[_kmer_length] = 0;
		fputs(str, out_file);

		string Kmer = string(str);
		if (Kmer <= PrevKmer)
			{
			Log("%s  prev\n", PrevKmer.c_str());
			Log("%s  next\n", Kmer.c_str());
			Die("Kmer <= PrevKmer");
			}
		PrevKmer = Kmer;

		if (counter <= UINT_MAX)
			fprintf(out_file, "\t%d\n", uint32(counter));
		else
			fprintf(out_file, "\t*\n");
		}

	CloseStdioFile(out_file);
	kmer_data_base.Close();
	}
