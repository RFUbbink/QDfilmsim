#pragma once
#include "Discontinuous.h"
#include <vector>
#include "Config.h"

class DisXcon:
	public DisCell
{
	
	enum contaminationType
	{
		contamination_X,
		contamination_Xmin,
		contamination_Xswitch
	};

private:
	array_type m_contaminations{};
	array_type m_contaminationCurrents{};

	const double m_currentConstantContaminantf{};
	const double m_currentConstantContaminants1{};
	const double m_currentConstantContaminants2{};

public:
	DisXcon(settings_array& settings);

	virtual void initializeConcentrations(double contaminantConcentration);
	virtual void calculatePotentialProfile();
	inline double negativeCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon);
	inline double negativeCurrente(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon);
	inline double positiveCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon);
	inline double neutralCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double eCon);
	void react();
	virtual void calculateCurrents();
	virtual void updateConcentrations();
	virtual void midSave(std::ofstream& midf);
};

