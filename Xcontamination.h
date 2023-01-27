#pragma once
#include "Electrochemistry.h"
#include <vector>
#include "Config.h"


class Xcontamination :
    public Cell
{
	enum contaminationType
	{
		contamination_X,
		contamination_Xplus,
		contamination_Xswitch
	};

private:
    array_type m_contaminations{};
    array_type m_contaminationCurrents{};

	const double m_currentConstantContaminant{};

public:
	Xcontamination(settings_array& settings);

	virtual void initializeConcentrations(double saltConcentration, double contaminantConcentration);
	virtual void calculatePotentialProfile();
	inline double negativeCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField);
	inline double negativeCurrente(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField);
	inline double positiveCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField);
	inline double neutralCurrent(const double concentrationLeft, const double concentrationRight, const double curCon);
	void react();
	virtual void calculateCurrents();
	virtual void updateConcentrations();
	virtual void midSave(std::ofstream& midf);
};

