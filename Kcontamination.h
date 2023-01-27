#pragma once
#include "Electrochemistry.h"
#include <vector>
#include "Config.h"

class Kcontamination :
    public Cell
{
	enum contaminationType
	{
		contamination_K,
		contamination_QDKmin,
		contamination_switch_maxCon
	};

private:
	array_type m_contaminations{};
	array_type m_contaminationCurrents{};

	const double m_currentConstantContaminant{};
	const double m_maximumElectronConcentration{};

public:
	Kcontamination(settings_array& settings);

	virtual void injectElectrons();
	virtual void initializeConcentrations(double saltConcentration, double contaminantConcentration);
	virtual void calculatePotentialProfile();
	inline double negativeCurrentmax(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double maximum = 1e27);
	inline double neutralCurrent(const double concentrationLeft, const double concentrationRight, const double curCon);
	void react();
	virtual void calculateCurrents();
	virtual void updateConcentrations();
	virtual void midSave(std::ofstream& midf);

};

