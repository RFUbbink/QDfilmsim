#pragma once
#include "Discontinuous.h"
#include <vector>
#include "Config.h"

class DisOxGM :
    public DisCell
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
	const double m_currentConstantXFilm{};
	const double m_currentConstantXSolution1{};
	const double m_currentConstantXSolution2{};
	const double m_k{};
	const double m_E0{};
	int m_stepCounter{};


public:
	DisOxGM(settings_array& settings, DOS_array& RR);

	inline double negativeCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon);
	inline double negativeCurrente(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon);
	inline double positiveCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon);
	inline double neutralCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double eCon);
	virtual void initializeConcentrations(double contaminantConcentration);
	virtual void calculatePotentialProfile();
	void updateRatesGM();
	void updateRatesBV();
	void reactGM();
	void reactBV();
	virtual void calculateCurrents();
	virtual void updateConcentrations();
	virtual void midSave(std::ofstream& midf);
};

