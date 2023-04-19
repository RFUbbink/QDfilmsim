#pragma once
#include <vector>
#include "Config.h"

class Cell
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
		carrier_electrons,
		carrier_cations,
		carrier_anions
	};

	array_type m_electrostatic{};
	array_type m_concentrations{};
	array_type m_currents{};

	double m_appliedBias{};									//appliedBias is the TOTAL potential drop over the entire system == the voltage difference between the working and counter electrode
	double m_currentCumulative{};							//Current counter
	double m_leakCurrentCumulative{}; 
	const double m_voltageIncrement{}; 
	double m_saltConcentration{};							//ion concentration in the electrolyte solution

	//space
	const array_type::size_type m_size{}; 					//the size of the cell array
	const array_type::size_type m_interfacePoint{};			//the position of the film/solution interface in the array
	const array_type::size_type m_referencePoint{};			//the position of the refernce electrode in the array
	const double m_referencePositionRelative{};				//reference postion/cellthickness, used in potential calculation
	const double m_thickness{};								//distance between the WE and Ce
	const double m_dx{};									//the size of a single cell

	//band
	const double m_injectionBarrier{};						//the injection barrier for electrons
	const double m_LUMO{};									//legacy constant I think, energy level of the CB edge
	double m_negativeElectrodeWF{};							//Fermi level in the WE material = work function of that material
	const double m_QDFillFactor{};							//the fraction of space that is filled by QD in the film, the rest is electrolyte
	const double m_electronEnergyFactor{};					//The parameter that allows to compensate the drif-diffusion currents for the fact that the more electrons you put in the CB, the more energy is required.

	//speedy constants
	const double m_currentConstantElectrons{};				//these are all used in the calculation of the currents, they are precalculated here to avoid unnecessary multiplications in the main code
	double m_currentConstantCationsFilm{};
	const double m_currentConstantCationsSolution{};
	double m_currentConstantAnionsFilm{};
	const double m_currentConstantAnionsSolution{};
	const double m_energyConvert{};
	const double m_energyConvertx{};
	const double m_poissonConstantFilm{};
	const double m_poissonConstantSolution{};
	const double m_currentConvert{};
	const double m_OhmicDropConstant{};
	double m_REcorrection{};
	int m_counter{};

public:

	Cell(settings_array& settings);

	double getCurrent();
	double getLeakCurrent();
	double getVoltage();

	virtual void injectElectrons(const DOS_array& DOS);
	virtual void calculatePotentialProfile();
	virtual void initializeConcentrations();
	inline double negativeCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField);
	inline double negativeCurrente(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField);
	inline double positiveCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField);
	virtual void calculateCurrents();
	virtual void updateConcentrations();
	void loadState();
	void resetInjection();
	void changeBias(double vbias);;
	Cell& operator++();
	Cell& operator--();
	virtual void midSave(std::ofstream& midf);
	virtual ~Cell() = default;
};

