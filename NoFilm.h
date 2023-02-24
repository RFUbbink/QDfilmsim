/*
This class is used to perform pure electrocal simulations on reductant/oxidant pairs. 
No quantum dot film is present. No discontinuous space intervals are employed.
*/


#pragma once
#include <vector>
#include "Config.h"
#include <iostream>

class NoFilmCell
{
protected: //all the constants here are calculated from the config upon initialization

	enum electrostaticType
	{
		es_secondDerivative,
		es_electricField,
		es_potential
	};

	enum carrierType
	{
		carrier_X,
		carrier_Xmin,
		carrier_cations,
		carrier_anions
	};

	array_type m_electrostatic{};
	array_type m_concentrations{};
	array_type m_currents{};

	double m_appliedBias{};
	double m_currentCumulative{};
	const double m_voltageIncrement{};
	double m_saltConcentration{};
	double m_Xconcentration{};

	//space
	static const array_type::size_type m_size{ precompiled::amountOfCells }; 	//the size of the cell array
	const array_type::size_type m_referencePoint{};			//the position of the refernce electrode in the array
	const double m_referencePositionRelative{};				//reference postion/cellthickness, used in potential calculation
	const double m_thickness{};
	const double m_dx{};									//the size of a single cell of the film

	//band
	const double m_E0{};

	//speedy constants
	const double m_currentConstantX{};				//these are all used in the calculation of the currents
	const double m_currentConstantCations{};
	const double m_currentConstantAnions{};
	const double m_energyConvert{};
	const double m_energyConvertx{};
	const double m_poissonConstantSolution{};
	const double m_currentConvert{};
	const double m_Nernst{};

public:

	NoFilmCell(settings_array& settings);

	double getCurrent();
	double getVoltage();
	double getLeakCurrent();

	void resetInjection();
	void injectElectrons(const DOS_array& DOS);
	void calculatePotentialProfile();
	void initializeConcentrations();
	inline double negativeCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon);
	inline double neutralCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double eCon);
	inline double positiveCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon);
	virtual void calculateCurrents();
	virtual void updateConcentrations();
	NoFilmCell& operator++();
	NoFilmCell& operator--();
	virtual void midSave(std::ofstream& midf);
	virtual ~NoFilmCell() = default;
};

