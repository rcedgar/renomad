#include "myutils.h"
#include "params.h"

/***
We proceed with p-value calculation only for:
	anchors with more than 50 total counts across all samples. 
		g_MinAnchorTotalBeforeFiltering = 50

We discard:
	anchors that have only one unique target
		g_MinTargetsPerAnchor = 2
		--anchor_unique_targets_threshold [1]

	anchors that appear in only 1 sample
		g_MinTargetsPerAnchor = 2
		--anchor_samples_threshold [1]

	(anchor, sample) pairs that have fewer than 6 counts.
		g_MinSampleTotal = 6
		--anchor_sample_counts_threshold [1]
	
Finally, we retain only anchors having more than 30 total 
counts after above thresholds were applied.
		// g_MinAnchorTotalAfterFiltering = // 30 NOT IMPLEMENTED
		--anchor_count_threshold [30]
***/

uint g_MinAnchorTotalBeforeFiltering = 50;
uint g_MinSamplesPerAnchor = 2;
uint g_MinTargetsPerAnchor = 2;
uint g_MinSampleTotal = 1; // says 6 in paper

// RCE filters
uint g_MinTargetTotal = 1;
uint g_MinTargetMaxCount = 1;

uint g_SplitCount = 50;
uint g_FuncCount = 10;
uint g_RetryCount = 1;
double g_MaxP = 0.05;
bool g_DoBonferroni = true;

uint g_InitialMaxTargetCount = 1000;

void SetParams_CmdLine()
	{
	g_MaxP = 0.5;
	if (optset_maxp)
		g_MaxP = opt(maxp);
	asserta(g_MaxP >= 0 && g_MaxP <= 1);

	g_RetryCount = 1;
	if (optset_retries)
		g_RetryCount = opt(retries);
	}

void SetParams_Struct(const NomadParams &NP)
	{
#define s(x)	g_##x = NP.x

	s(MinAnchorTotalBeforeFiltering);
	s(MinSamplesPerAnchor);
	s(MinTargetsPerAnchor);
	s(MinSampleTotal);
	s(MinTargetMaxCount);
	s(SplitCount);
	s(FuncCount);
	s(RetryCount);
	s(DoBonferroni);

#undef s
	}
