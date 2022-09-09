#pragma once

extern uint g_MinAnchorTotalBeforeFiltering;
extern uint g_MinTargetsPerAnchor;
extern uint g_MinSamplesPerAnchor;
extern uint g_MinSampleTotal;

// RCE filters
extern uint g_MinTargetTotal;
extern uint g_MinTargetMaxCount;

extern uint g_SplitCount;
extern uint g_FuncCount;
extern uint g_RetryCount;
extern double g_MaxP;

extern bool g_DoBonferroni;

extern uint g_InitialMaxTargetCount;

struct NomadParams
	{
	uint MinAnchorTotalBeforeFiltering = 50;
	uint MinSamplesPerAnchor = 2;
	uint MinTargetsPerAnchor = 2;
	uint MinSampleTotal = 1;
	uint MinTargetMaxCount = 1;
	uint SplitCount = 50;
	uint FuncCount = 10;
	uint RetryCount = 1;
	bool DoBonferroni = true;
	};

void SetParams_CmdLine();
void SetParams_Struct(const NomadParams &NP);
