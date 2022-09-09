#pragma once

#include "atab.h"

const uint LONGVEC_LENGTH = 1000003;

class Nomad : public Atab
	{
public:
	uint m_SampleCount = 0;
	uint m_MaxTargetCount = 0;

	uint *m_SelectedSamples = 0;
	uint *m_SelectedTargets = 0;

	uint m_SelectedSampleCount = 0;
	uint m_SelectedTargetCount = 0;

// Per-target, per-iteration fs
//   represented as long vector
	bool *m_f_vec = 0;

// Per-sample total count
	uint *m_njs = 0;
	double *m_sqrt_njs = 0;

// Per-sample count w.r.t. fs
	double *m_mujs = 0;

// Per-sample S
	double *m_Sjs = 0;

// Per-sample, per-iteration partition +/-1
//   represented as long vector
	int *m_c_vec = 0;

// Total of counts in matrix
	uint m_M = 0;

// Total of function-selected counts
	double m_m = 0;
	double m_mu = 0;

	bool *m_fs = 0;
	int *m_cs = 0;

public:
	void Clear();

	void Init(uint ka, uint kt, uint SampleCount,
	  uint InitialMaxTargetCount);

	void Init(uint ka, uint kt, 
	  const vector<string> &SampleNames,
	  uint InitialMaxTargetCount);

	bool Select();

	void Set_fs();
	void Set_cs();
	void Set_njs();
	void Set_mujs();
	void Set_Sjs();
	double CalcP1();
	double CalcP();
	};
