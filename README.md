# RENOMAD
Implementation of NOMAD (https://doi.org/10.1101/2022.06.24.497555)

### K-mer counting with KMC3

First run KMC3 for each sample, here is an example which assumes `ERR5296631_1.fastq` is  present in the current directory.

<pre>
SRA=ERR5296631
$kmc \
  -k54 \
  -m32 \
  -b \
  -jkmc_json/${SRA}.json \
  fq/${SRA}_1.fastq \
  kmc_out/${SRA}.kmc \
  kmc_tmp/
</pre>

### Sort KMC3 output files into lexicographic order.

<pre>
SRA=ERR5296631
kmc_tools transform $SRA.kmc sort $SRA.sorted.kmc
</pre>

### Run renomad

Create a text file named `kmcfiles.txt` with pathnames of the sorted KMC output files, one per line. Then run `renomad` as follows:

<pre>
renomad \
  -joinp files.txt \
  -maxp 0.05 \
  -tsv3out output.tsv
</pre>

The `maxp` option is the maximum P value to report an anchor. Default is 0.05.

Output format is tab-separated text with three fields: (1) count, (2) anchor+target sequence, (3) sample_name. The sample name is extracted from the KMC filename by stripping the path name and extension. 
