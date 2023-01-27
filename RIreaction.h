#pragma once
#include "Electrochemistry.h"
class RIreaction :
    public Cell
{

	enum contaminants
	{
		contaminant_X,
		contaminant_Xmin,
		switch_X
	};


private:
	array_type m_contaminants{};
	array_type m_contaminantCurrents{};

	const double m_contaminantConcentration{};
	const double m_currentConstantX{};
	const double m_k1forward{};
	const double m_k1backward{};
	const double m_k2{};
	const double m_maxTrapConcentration{};


public:
	RIreaction(settings_array& settings);

	virtual void initializeConcentrations(double contaminantConcentration);
	virtual void injectElectrons(const DOS_array& DOS);
	virtual void calculatePotentialProfile();
	void react();
	inline double negativeCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField);
	inline double negativeCurrente(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField);
	inline double positiveCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField);
	double neutralCurrent(const double concentrationLeft, const double concentrationRight, const double curCon);
	virtual void calculateCurrents();
	virtual void updateConcentrations();
	virtual void midSave(std::ofstream& midf);
};

