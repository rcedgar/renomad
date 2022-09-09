#include "myutils.h"
#include "nomad.h"
#include "params.h"

#define TRACE	0

void Nomad::Clear()
	{
	Atab::Clear();

	freez(m_SelectedSamples);
	freez(m_SelectedTargets);
	freez(m_f_vec);
	freez(m_njs);
	freez(m_mujs);
	freez(m_Sjs);
	freez(m_c_vec);

	m_SampleCount = 0;
	m_MaxTargetCount = 0;
	m_SelectedSampleCount = 0;
	m_SelectedTargetCount = 0;
	m_M = 0;
	m_m = 0;
	m_fs = 0;
	m_cs = 0;
	}

void Nomad::Init(uint ka, uint kt, uint SampleCount,
  uint InitialMaxTargetCount)
	{
	vector<string> SampleNames;
	for (uint i = 0; i < SampleCount; ++i)
		{
		string SampleName;
		Ps(SampleName, "Sam%u", i+1);
		SampleNames.push_back(SampleName);
		}
	Init(ka, kt, SampleNames, InitialMaxTargetCount);
	}

void Nomad::Init(uint ka, uint kt,
  const vector<string> &SampleNames,
  uint InitialMaxTargetCount)
	{
	Clear();
	Atab::Init(ka, kt, SampleNames, InitialMaxTargetCount);

#if TRACE
	Log("\n");
	Log("Nomad::Init(ka=%u, kt=%u, IMTC=%u)\n",
	  ka, kt, InitialMaxTargetCount);
#endif
	m_SampleCount = SIZE(SampleNames);

	m_SelectedSamples = myalloc(uint, m_SampleCount);
	m_SelectedTargets = myalloc(uint, InitialMaxTargetCount);

	m_njs = myalloc(uint, m_SampleCount);
	m_sqrt_njs = myalloc(double, m_SampleCount);
	m_mujs = myalloc(double, m_SampleCount);
	m_Sjs = myalloc(double, m_SampleCount);;

	m_f_vec = myalloc(bool, LONGVEC_LENGTH);
	m_c_vec = myalloc(int, LONGVEC_LENGTH);

	for (uint i = 0; i < LONGVEC_LENGTH; ++i)
		{
		uint r = randu32()%2;
		m_f_vec[i] = (r == 1);
		}

	for (uint i = 0; i < LONGVEC_LENGTH; ++i)
		{
		uint r = randu32()%2;
		m_c_vec[i] = (r == 1 ? 1 : -1);
		}
	}

void Nomad::Set_fs()
	{
	const uint TargetCount = m_TargetCount;
	uint N = (LONGVEC_LENGTH - TargetCount);
	uint Offset = randu32()%N;
	m_fs = m_f_vec + Offset;
#if TRACE
	{
	Log("Nomad::Set_fs() ");
	for (uint ii = 0; ii < m_SelectedTargetCount; ++ii)
		{
		uint TargetIndex = m_SelectedTargets[ii];
		Log("%c", tof(m_fs[TargetIndex]));
		}
	Log("\n");
	}
#endif
	}

void Nomad::Set_cs()
	{
	uint N = (LONGVEC_LENGTH - m_SampleCount);
	uint Offset = randu32()%N;
	m_cs = m_c_vec + Offset;
#if TRACE
	{
	Log("Nomad::Set_cs()");
	for (uint jj = 0; jj < m_SelectedSampleCount; ++jj)
		{
		uint SampleIndex = m_SelectedSamples[jj];
		Log(" %+d", m_cs[SampleIndex]);
		}
	Log("\n");
	}
#endif
	}

void Nomad::Set_njs()
	{
	m_M = 0;
	const uint TargetCount = m_TargetCount;
	for (uint jj = 0; jj < m_SelectedSampleCount; ++jj)
		{
		const uint SampleIndex = m_SelectedSamples[jj];
		uint nj = 0;
		for (uint ii = 0; ii < m_SelectedTargetCount; ++ii)
			{
			const uint TargetIndex = m_SelectedTargets[ii];
			nj += m_Counts[TargetIndex][SampleIndex];
			}

		m_njs[SampleIndex] = nj;
		m_sqrt_njs[SampleIndex] = sqrt(nj);
		m_M += nj;
		}
#if TRACE
	{
	Log("Nomad::Set_njs()");
	for (uint jj = 0; jj < m_SelectedSampleCount; ++jj)
		{
		uint SampleIndex = m_SelectedSamples[jj];
		Log(" %u", m_njs[SampleIndex]);
		}
	Log("\n");
	}
#endif
	}

void Nomad::Set_mujs()
	{
	m_m = 0;
	const uint TargetCount = m_TargetCount;
	for (uint jj = 0; jj < m_SelectedSampleCount; ++jj)
		{
		uint SampleIndex = m_SelectedSamples[jj];
		uint m_j = 0;

		for (uint kk = 0; kk < m_SelectedTargetCount; ++kk)
			{
			uint TargetIndex = m_SelectedTargets[kk];
			if (m_fs[TargetIndex])
				m_j += m_Counts[TargetIndex][SampleIndex];
			}

		double mu_j = 0;
		if (m_j > 0)
			{
			uint nj = m_njs[SampleIndex];
			mu_j = double(m_j)/nj;
			m_m += m_j;
			}
		m_mujs[SampleIndex] = mu_j;
		}
	asserta(m_m <= m_M);
	m_mu = 0;
	if (m_M > 0)
		{
		m_mu = double(m_m)/double(m_M);
		asserta(m_mu <= 1);
		}
#if TRACE
	{
	Log("Nomad::Set_mujs()");
	for (uint jj = 0; jj < m_SelectedSampleCount; ++jj)
		{
		uint SampleIndex = m_SelectedSamples[jj];
		Log(" %.2g", m_mujs[SampleIndex]);
		}
	Log("\n");
	Log("mu = %.2g\n", m_mu);
	}
#endif
	}

void Nomad::Set_Sjs()
	{
	for (uint jj = 0; jj < m_SelectedSampleCount; ++jj)
		{
		uint SampleIndex = m_SelectedSamples[jj];
		double sqrt_n_j = m_sqrt_njs[SampleIndex];
		double mu_j = m_mujs[SampleIndex];
		double S_j = sqrt_n_j*(mu_j - m_mu);
		m_Sjs[SampleIndex] = S_j;
		}
#if TRACE
	{
	Log("Nomad::Set_Sjs()");
	for (uint jj = 0; jj < m_SelectedSampleCount; ++jj)
		{
		uint SampleIndex = m_SelectedSamples[jj];
		Log(" %.3g", m_Sjs[SampleIndex]);
		}
	Log("\n");
	}
#endif
	}

bool Nomad::Select()
	{
	m_SelectedSampleCount = 0;
	m_SelectedTargetCount = 0;

	if (m_Total < g_MinAnchorTotalBeforeFiltering)
		return false;
	if (m_TargetCount < g_MinTargetsPerAnchor)
		return false;

	for (uint SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
		{
		uint SampleTotal = m_SampleTotals[SampleIndex];
		if (SampleTotal >= g_MinSampleTotal)
			{
			m_SelectedSamples[m_SelectedSampleCount] = SampleIndex;
			++m_SelectedSampleCount;
			}
		}
	if (m_SelectedSampleCount < g_MinSamplesPerAnchor)
		return false;

	for (uint TargetIndex = 0; TargetIndex < m_TargetCount; ++TargetIndex)
		{
		uint TargetTotal = m_TargetTotals[TargetIndex];
		if (TargetTotal >= g_MinTargetTotal)
			{
			if (m_SelectedTargetCount >= g_InitialMaxTargetCount)
				break;
			m_SelectedTargets[m_SelectedTargetCount] = TargetIndex;
			++m_SelectedTargetCount;
			}
		}
	if (m_SelectedTargetCount < g_MinTargetsPerAnchor)
		return false;

#if TRACE
	Log("Nomad::Select()\n");
	Log("  Samples: ");
	for (uint jj = 0; jj < m_SelectedSampleCount; ++jj)
		Log(" %d", m_SelectedSamples[jj]);
	Log("\n");
	Log("  Targets: ");
	for (uint ii = 0; ii < m_SelectedTargetCount; ++ii)
		Log(" %d", m_SelectedTargets[ii]);
	Log("\n");
#endif
	return true;
	}

double Nomad::CalcP1()
	{
	double S = 0;
	double Sum_cj2_ngt0 = 0;
	double Sum_cj_sqrtnj = 0;

	for (uint jj = 0; jj < m_SelectedSampleCount; ++jj)
		{
		uint SampleIndex = m_SelectedSamples[jj];
		uint n_j = m_njs[SampleIndex];
		int c_j = m_cs[SampleIndex];
		double S_j = m_Sjs[SampleIndex];
		double sqrt_nj = m_sqrt_njs[SampleIndex];

		S += S_j*c_j;
		Sum_cj_sqrtnj += c_j*sqrt_nj;
		if (n_j > 0)
			Sum_cj2_ngt0 += 1; // c_j*c_j;
		}
	const double Sum_cj2 = m_SampleCount;
	if (Sum_cj_sqrtnj == 0)
		return 1.0;

	double ablob2 = float(m_M)*Sum_cj2/float(Sum_cj_sqrtnj*Sum_cj_sqrtnj);
	double ablob = sqrt(ablob2);
	double a = 1.0/(1 + ablob);
	double exp1_top = 2.0*(1 - a*a)*S*S;
	double exp1 = exp1_top/Sum_cj2_ngt0;
	double exp2_top = 2*a*a*m_M*S*S;
	double exp2 = exp2_top/(Sum_cj_sqrtnj*Sum_cj_sqrtnj);

	double P = 2*exp(-exp1) + 2*exp(-exp2);
	asserta(!isnan(P));
	if (P > 1.0)
		P = 1.0;
	return P;
	}

double Nomad::CalcP()
	{
	bool Ok = Select();
	if (!Ok)
		return 1.0;

	Set_njs();

	uint TestCount = 0;
	double MinP = 1.0;
	double LastP = 1.0;
	for (uint Try = 0; Try < g_RetryCount; ++Try)
		{
		for (uint SplitIndex = 0; SplitIndex < g_SplitCount; ++SplitIndex)
			{
			++TestCount;
			Set_cs();
			for (uint FuncIndex = 0; FuncIndex < g_FuncCount; ++FuncIndex)
				{
#if TRACE
				Log("\n");
				Log("=== Split %u/%u, func %u/%u ===\n",
				  SplitIndex, g_SplitCount, FuncIndex, g_FuncCount);
#endif
				Set_fs();
				Set_mujs();
				Set_Sjs();
				double P1 = CalcP1();
				if (P1 < MinP)
					MinP = P1;
#if TRACE
				Log("P1 = %.3g, MinP = %.3g\n", P1, MinP);
#endif
				}
			}
		double P = MinP*TestCount;
		if (P <= 0.05)
			return P;
		if (P >= LastP)
			return 1.0;
		}

	double P = MinP*TestCount;
	if (P > 1.0)
		P = 1.0;
	return P;
	}
