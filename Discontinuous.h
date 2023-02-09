/*
* This is now the most important class, which does everything except the electrochemical side reactions part. DisOxGM does that and can also be used to simulate normal curves but is slower than this one due to more steps and more memory.
* 
*/

#pragma once
#include <vector>
#include "Config.h"
#include <iostream>

class DisCell
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

	double m_appliedBias{};									//the potential vs the RE
	double m_currentCumulative{};							//The amount of electrons that entered the film since the last time it was reported
	double m_leakCurrentCumulative{};						//same for the leak electrons, not used that much
	const double m_voltageIncrement{};						//amount the voltage is incremented every voltage step
	double m_saltConcentration{};							//concentration of cation and anions in the electrolyte solution at the start. 

	//space
	static const array_type::size_type m_size{ precompiled::amountOfCells }; 	//the size of the cell array
	static const array_type::size_type m_Isize{ precompiled::InterfaceCells }; 	//the size of the interface
	const array_type::size_type m_interfacePoint{};			//the position of the film/solution interface in the array
	const array_type::size_type m_referencePoint{};			//the position of the refernce electrode in the array
	const double m_referencePositionRelative{};				//reference postion/cellthickness, used in potential calculation
	const double m_thickness{};
	const double m_dxf{};									//the size of a single cell of the film
	const double m_dxs1{};									//the size of a single cell of the solution
	const double m_dxs2{};									//the size of a single cell of the solution

	//band
	const double m_injectionBarrier{};						//the injection barrier for electrons
	const double m_densityOfStates{};						//the effective density of states for electrons in the QD film
	const double m_LUMO{};									//LUMO of the QDs (not used as it is also in the DOS array)
	double m_negativeElectrodeWF{};							//work function of the working electrode
	const double m_QDFillFactor{};							//the fraction of space that is filled by QD in the film, the rest is electrolyte
	const double m_electronEnergyFactor{};					//The secret electron energy factor which is also in the DOS array, corrects the drift/diffusion current of electrons for the the difference in energy state.

	//speedy constants
	const double m_currentConstantElectrons{};				//these are all used in the calculation of the currents, look in the constructor if you want to know how they are calculated. 
	const double m_currentConstantElectrons2{};				//Used to vaoid calculations inside the functions
	double m_currentConstantCationsFilm{};
	double m_currentConstantCationsFilm2{};
	const double m_currentConstantCationsSolution1{};
	const double m_currentConstantCationsSolution2{};
	double m_currentConstantAnionsFilm{};
	double m_currentConstantAnionsFilm2{};
	const double m_currentConstantAnionsSolution1{};
	const double m_currentConstantAnionsSolution2{};
	const double m_energyConvert{};
	const double m_energyConvertxf{};
	const double m_energyConvertxs1{};
	const double m_energyConvertxs2{};
	const double m_poissonConstantFilm{};
	const double m_poissonConstantSolution{};
	const double m_currentConvert{};

public:

	DisCell(settings_array& settings);

	double getCurrent();
	double getLeakCurrent();
	double getVoltage();

	virtual void injectElectrons(const DOS_array& DOS);
	virtual void calculatePotentialProfile();
	virtual void initializeConcentrations();
	inline double negativeCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon);
	inline double negativeCurrente(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon);
	inline double positiveCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon);
	virtual void calculateCurrents();
	virtual void updateConcentrations();
	void changeBias(double vbias);
	void resetInjection();
	void adjustMobility(settings_array& settings);
	DisCell& operator++();
	DisCell& operator--();
	virtual void midSave(std::ofstream& midf);
	virtual ~DisCell() = default;
};

