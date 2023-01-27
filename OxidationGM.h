#pragma once
#include "Electrochemistry.h"
#include <vector>
#include "Config.h"

class OxidationGM :
    public Cell
{
	enum contaminants
	{
		contaminant_X,
		contaminant_Xmin,
		switch_X
	};

	enum rates
	{
		rate_forward,
		rate_backward,
		rate_NONVALID
	};

private:
	array_type m_contaminants{};
	array_type m_contaminantCurrents{};
	array_type m_reactionRates{};

	DOS_array m_RRlookup{};

	const double m_contaminantConcentration{};
	const double m_currentConstantX{};
	const double m_k{};
	int m_stepCounter{};


public:
	OxidationGM(settings_array& settings, DOS_array& RR);

	virtual void initializeConcentrations(double contaminantConcentration);
	virtual void calculatePotentialProfile();
	void updateRates();
	void react();
	inline double negativeCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField);
	inline double negativeCurrente(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField);
	inline double positiveCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField);
	double neutralCurrent(const double concentrationLeft, const double concentrationRight, const double curCon);
	virtual void calculateCurrents();
	virtual void updateConcentrations();
	virtual void midSave(std::ofstream& midf);
};
