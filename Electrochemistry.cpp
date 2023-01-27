//The most basic cell object. Has continuous lamella spacing.

#include <cmath>
#include <cstddef>
#include <vector>
#include <array>
#include <numeric>
#include <iostream>
#include <fstream>
#include "Config.h"
#include "Electrochemistry.h" 

using array_type = std::vector<std::vector<double>>;

//Member functions of Cell class:

Cell::Cell(settings_array& settings) //initialize all the constants necessary for operation, using settings array
//an estimation of the initial voltage drop over the counterelectrode is made to calculate the initial appliedBias.
	:m_appliedBias{ settings[s_startVoltage] + settings[s_startVoltage] * settings[s_workingElectrodeArea]
						/ settings[s_counterElectrodeArea] * settings[s_epsilonrFilm] / settings[s_epsilonrSolution] },
	m_voltageIncrement{ settings[s_voltageIncrement] },
	m_maximumIonConcentration{ settings[s_maximumIonConcentration] },
	m_saltConcentration{ settings[s_ionConcentration] },
	m_interfacePoint{ static_cast<array_type::size_type>(1 + static_cast<int>(m_size * settings[s_filmThickness] / settings[s_cellThickness])) },
	m_referencePoint{ static_cast<array_type::size_type>(m_size * settings[s_refPosition] / settings[s_cellThickness]) },
	m_referencePositionRelative{ settings[s_refPosition] / settings[s_cellThickness] },
	m_thickness{ settings[s_cellThickness] },
	m_dx{ settings[s_cellThickness] / (m_size - 1) },
	m_injectionBarrier{ settings[s_LUMO] - settings[s_negativeElectrodeWF] },
	m_densityOfStates{ settings[s_densityOfStates] },
	m_LUMO{ settings[s_LUMO] },
	m_negativeElectrodeWF{ settings[s_negativeElectrodeWF] },
	m_QDFillFactor{ 1 - settings[s_QDspacefill] },
	m_electronEnergyFactor{ settings[39] / (settings[s_cellThickness] / (m_size - 1)) },
	m_currentConstantElectrons{ settings[s_electronMobility] * settings[s_dt] / settings[s_cellThickness] * m_size },
	m_currentConstantCationsFilm{ settings[s_cationMobilityFilm] * settings[s_dt] / settings[s_cellThickness] * m_size },
	m_currentConstantCationsSolution{ settings[s_cationMobilitySolution] * settings[s_dt] / settings[s_cellThickness] * m_size },
	m_currentConstantAnionsFilm{ settings[s_anionMobilityFilm] * settings[s_dt] / settings[s_cellThickness] * m_size },
	m_currentConstantAnionsSolution{ settings[s_anionMobilitySolution] * settings[s_dt] / settings[s_cellThickness] * m_size },
	m_energyConvert{ phys::k * settings[s_temperature] / phys::q },
	m_energyConvertx{ phys::k * settings[s_temperature] / phys::q / settings[s_cellThickness] * m_size },
	m_poissonConstantFilm{ -phys::q / phys::eps0 / settings[s_epsilonrFilm] },
	m_poissonConstantSolution{ -phys::q / phys::eps0 / settings[s_epsilonrSolution] },
	m_currentConvert{ phys::q * settings[s_cellThickness] / m_size / settings[s_dt] / 10000 },
	m_OhmicDropConstant{ 2e-5 * settings[s_cellThickness] / m_size / settings[s_dt] / (settings[s_ionConcentration] * (settings[s_anionMobilitySolution] + settings[s_cationMobilitySolution]) / 2) }, // Ohmic drop distance = 0.2 mm
	m_REcorrection{ 0 },
	m_counter{0}
{
	m_electrostatic = array_type(3, std::vector<double>(m_size));		//create the potential array
	m_electrostatic[es_potential][0] = settings[s_startVoltage];					//apply the inital bias: updated everytime there is a voltage increment.
																		//This array value should be constant unless there is an increment
	m_concentrations = array_type(3, std::vector<double>(m_size));		//create the concentration array
	m_currents = array_type(3, std::vector<double>(m_size));			//create the current array
}

double Cell::getCurrent() //returns amount of electrons that have entered the film from the working electrode since the last time this function was called
{
	double curcur{ m_currentCumulative * m_currentConvert };
	m_currentCumulative = 0;
	return curcur;
}

double Cell::getLeakCurrent() //similar to getCurrent but returns the amount of electrons that have eneterd the counterelectrode
{
	double curcur{ m_leakCurrentCumulative * m_currentConvert };
	m_leakCurrentCumulative = 0;
	return curcur;
}

double Cell::getVoltage() //Returns the current applied potential
{
	return m_electrostatic[es_potential][0];
}

void Cell::injectElectrons(const DOS_array& DOS)
{
	//These lines use the real DOS and Fermi-Dirac distribution to calculate the electron concentration.
	//DOS array also holds the "official" starting energy (DOS[0]) and energy integrating step dE (DOS[1])
	double n{};
	double energy{ DOS[0]};
	double dE{ DOS[1]};
	//integrating from Estart to the end of DOS eV using steps of 0.01 eV
	for (DOS_array::size_type i{ 3 }; i < 300; ++i)
	{
		n += DOS[i] / (exp((energy - (m_negativeElectrodeWF - m_electrostatic[es_potential][0]+ m_electrostatic[es_potential][1])) / m_energyConvert) + 1)* dE;
		energy += dE;
	}
	
	//The line below adjusts the electron concentration in the working electrode up or down
	//Over time, this leads to equilibrium where the concentration in the first cell of the film is equal to the one just calculated using the DOS
	//So, iteratively, chemical equilibrium is enforced..
	m_concentrations[carrier_electrons][0] = m_concentrations[carrier_electrons][0]*n/ m_concentrations[carrier_electrons][1];
}

void Cell::calculatePotentialProfile() //This is the most critical function in the sim. It calculates the electrostatic potential by solving the Poisson equation. I could not use a normal ODE I think, but I may be wrong.
//I had to adapt the ODE so that it can apply the boundary conditions properly, and ended up doing this iterative approach. It typically can solve it in 4-6 iterations, but it may be more in the middle of the sim. 
{
	//calculate the second space derivative of the electrostatic potential for every cell based on the Poisson equation
	
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantFilm * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i] - m_concentrations[carrier_anions][i]);
	}

	for (array_type::size_type i{ m_interfacePoint }; i < (m_size-1); ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantSolution * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i] - m_concentrations[carrier_anions][i]);
	}

	double potAtReference{ 1.0 };
	while (std::abs(potAtReference) > 0.0001) //typically takes 2 rounds, very good!
	{
		//electricFieldStart is the initial guess for a boundary condition (assumes that the potential drops linearly over the cell at the start) because the real boundary conditions are impossible to apply directly 
		//This boundary condition will then be updated based on the potential profile that is calculated
		//This process is repeated untal a potential profile is generated that satisfies the real boundary conditions (pot at negative electrode == applied potential && pot at reference electrode == 0)
		double electricFieldStart{ (m_appliedBias) / m_thickness };
		double differenceEF{ 1.0 };
		double electricFieldIntegral{};
		while (std::abs(differenceEF) > 0.00000001) //Speed depends a bit, ranging from 2-5 rounds
		{
			m_electrostatic[es_electricField][0] = electricFieldStart; //apply guess boundary condition
			for (array_type::size_type i{ 0 }; i < (m_size - 2); ++i) //build electric field sequentially using said boundary condition
			{
				m_electrostatic[es_electricField][i + 1] = m_electrostatic[es_electricField][i] - m_electrostatic[es_secondDerivative][i + 1] * m_dx;
			}

			//check to see if real boundary condition is met (pot at negative elctrode == applied potential)
			electricFieldIntegral = std::accumulate(m_electrostatic[es_electricField].begin(), m_electrostatic[es_electricField].end(), 0.0) * m_dx;
			differenceEF = (m_appliedBias - electricFieldIntegral) / (m_appliedBias);
			//and apply the updated boundary condition for next round
			electricFieldStart = electricFieldStart * (differenceEF + 1);
		}

		for (array_type::size_type i{ 0 }; i < (m_size - 1); ++i) // update potential profile using the aquired electrid field
		{
			m_electrostatic[es_potential][i + 1] = m_electrostatic[es_potential][i] - m_electrostatic[es_electricField][i] * m_dx;
		}

		//check to see if real boundary condition is met (pot at reference electrode == 0)
		potAtReference = m_electrostatic[es_potential][m_referencePoint];
		//and update the bias to reflect this for the next round
		//divide by relative position of the reference point (default 1.5) for fastest conversion
		m_appliedBias += (potAtReference)/m_referencePositionRelative;
	}
}

void Cell::initializeConcentrations(double contaminantConcentration)
{
	//fills the concentration array with cations and anions
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_concentrations[carrier_cations][i] = m_saltConcentration* m_QDFillFactor;
		m_concentrations[carrier_anions][i] = m_saltConcentration* m_QDFillFactor;
	}

	for (array_type::size_type i{ m_interfacePoint }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_cations][i] = m_saltConcentration;
		m_concentrations[carrier_anions][i] =  m_saltConcentration;
	}
	contaminantConcentration;

	//Set the concentration of electrons to 1 (does not work for the injection function)
	m_concentrations[carrier_electrons][0] = 1;
	m_concentrations[carrier_electrons][1] = 1;
}

inline double Cell::injectionCurrent(const double concentrationElectrode, const double concentrationFilm)//legacy function
{
	return 0.001 * (concentrationElectrode - concentrationFilm);
}

inline double Cell::negativeCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (-concentrationLeft * electricField + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
}

inline double Cell::negativeCurrente(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	//if (concentrationLeft > 1e25)
	return (-concentrationLeft * (electricField - m_electronEnergyFactor*(concentrationLeft - concentrationRight)/(concentrationLeft+1)) + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
	//else
	//	return (-concentrationLeft * electricField  + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
}

inline double Cell::negativeCurrentmax(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double maximum)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	double current{ (-concentrationLeft * electricField + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon };
	return ((current < maximum) ? current : maximum);
}

inline double Cell::positiveCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField)
{
	//returns the amount of positively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (concentrationRight * electricField + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
}

inline double Cell::positiveCurrentmax(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double maximum)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	double current{ (concentrationRight * electricField + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon };
	return ((-current < maximum) ? current : -maximum);
}

void Cell::calculateCurrents()
{
	//electrons first, only up to the interface, using the current function to fill the current array
	//the first step is the injection current

	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_currents[carrier_electrons][i] = negativeCurrente(m_concentrations[carrier_electrons][i - 1], m_concentrations[carrier_electrons][i],
			m_currentConstantElectrons, m_electrostatic[es_electricField][i - 1]);
	}

	//cations next, needs two different loops for the film and solution sections. Ions cannot enter the electrodes, so we dont need to consider the first and last cell
	for (array_type::size_type i{ 2 }; i < m_interfacePoint; ++i)
	{
		m_currents[carrier_cations][i] = positiveCurrent(m_concentrations[carrier_cations][i - 1], m_concentrations[carrier_cations][i],
			m_currentConstantCationsFilm, m_electrostatic[es_electricField][i - 1]);
	}

	for (array_type::size_type i{ m_interfacePoint + 1 }; i < (m_size -1); ++i)
	{
		m_currents[carrier_cations][i] = positiveCurrent(m_concentrations[carrier_cations][i - 1], m_concentrations[carrier_cations][i],
			m_currentConstantCationsSolution, m_electrostatic[es_electricField][i - 1]);
	}

	//the anions which are very similar to cations
	for (array_type::size_type i{ 2 }; i < m_interfacePoint; ++i)
	{
		m_currents[carrier_anions][i] = negativeCurrent(m_concentrations[carrier_anions][i - 1], m_concentrations[carrier_anions][i],
			m_currentConstantAnionsFilm, m_electrostatic[es_electricField][i - 1]);
	}

	for (array_type::size_type i{ m_interfacePoint + 1 }; i < (m_size - 1); ++i)
	{
		m_currents[carrier_anions][i] = negativeCurrent(m_concentrations[carrier_anions][i - 1], m_concentrations[carrier_anions][i],
			m_currentConstantAnionsSolution, m_electrostatic[es_electricField][i - 1]);
	}
	
	//finally, we need to update the interface currents which we consider the current OVER the interface to be "slow"
	//this means limited by the lower concentration in the film, so we use the solution concentration*QDfillfactor to reflect this
	
	//LN(C) implementation
	//double catcur{ (sqrt(m_concentrations[carrier_cations][m_interfacePoint - 1] * m_concentrations[carrier_cations][m_interfacePoint] * m_QDFillFactor) * m_electrostatic[es_electricField][m_interfacePoint - 1] +
	//	m_energyConvertx * (m_concentrations[carrier_cations][m_interfacePoint - 1] - m_concentrations[carrier_cations][m_interfacePoint] * m_QDFillFactor)) * m_currentConstantCationsFilm };
	//m_currents[carrier_cations][m_interfacePoint] = (catcur < m_concentrations[carrier_cations][m_interfacePoint] ? catcur : m_concentrations[carrier_cations][m_interfacePoint]);
	//if (m_currents[carrier_cations][m_interfacePoint] < (m_concentrations[carrier_cations][m_interfacePoint - 1] - (m_concentrations[carrier_electrons][1] + 5e25)))
	//	m_currents[carrier_cations][m_interfacePoint] = (m_concentrations[carrier_cations][m_interfacePoint - 1] - (m_concentrations[carrier_electrons][1] + 5e25));
	
	//LN(A) implementation
	//double ancur{ (-sqrt(m_concentrations[carrier_anions][m_interfacePoint - 1] * m_concentrations[carrier_anions][m_interfacePoint]) * m_electrostatic[es_electricField][m_interfacePoint - 1] +
	//	m_energyConvertx * (m_concentrations[carrier_anions][m_interfacePoint - 1] - m_concentrations[carrier_anions][m_interfacePoint] * m_QDFillFactor)) * m_currentConstantAnionsFilm };
	//m_currents[carrier_anions][m_interfacePoint] = (ancur < m_concentrations[carrier_anions][m_interfacePoint] ? ancur : m_concentrations[carrier_anions][m_interfacePoint]);
	
	//double catcur{ (m_concentrations[carrier_cations][m_interfacePoint] * m_QDFillFactor * 3 * m_electrostatic[es_electricField][m_interfacePoint - 1] +
	//	m_energyConvertx * (m_concentrations[carrier_cations][m_interfacePoint - 1] - m_concentrations[carrier_cations][m_interfacePoint] * m_QDFillFactor)) * m_currentConstantCationsFilm };
	//m_currents[carrier_cations][m_interfacePoint] = (catcur < m_concentrations[carrier_cations][m_interfacePoint] ? catcur : m_concentrations[carrier_cations][m_interfacePoint]);

	//BOLTZMANN implementation
	//double catcur{ (m_concentrations[carrier_cations][m_interfacePoint - 1] - m_concentrations[carrier_cations][m_interfacePoint] * m_QDFillFactor *
	//												exp((m_electrostatic[es_potential][m_interfacePoint] - m_electrostatic[es_potential][m_interfacePoint - 1]) / m_energyConvert))/100000 };
	//m_currents[carrier_cations][m_interfacePoint] = (catcur < m_concentrations[carrier_cations][m_interfacePoint] ? catcur : m_concentrations[carrier_cations][m_interfacePoint]);
	//double ancur{(m_concentrations[carrier_anions][m_interfacePoint - 1] - m_concentrations[carrier_anions][m_interfacePoint] * m_QDFillFactor *
	//												exp((m_electrostatic[es_potential][m_interfacePoint - 1] - m_electrostatic[es_potential][m_interfacePoint]) / m_energyConvert))/100000 };
	//m_currents[carrier_anions][m_interfacePoint] = (ancur < m_concentrations[carrier_anions][m_interfacePoint] ? ancur : m_concentrations[carrier_anions][m_interfacePoint]);

	//REGULAR DD implementation
	m_currents[carrier_cations][m_interfacePoint] = positiveCurrent(m_concentrations[carrier_cations][m_interfacePoint - 1], m_concentrations[carrier_cations][m_interfacePoint] * m_QDFillFactor,
		m_currentConstantCationsFilm, m_electrostatic[es_electricField][m_interfacePoint - 1]);
	m_currents[carrier_anions][m_interfacePoint] = negativeCurrent(m_concentrations[carrier_anions][m_interfacePoint - 1], m_concentrations[carrier_anions][m_interfacePoint]*m_QDFillFactor,
		m_currentConstantAnionsFilm, m_electrostatic[es_electricField][m_interfacePoint - 1]);


	//update total electrons that have entered since the last getCurrent() call:
	m_currentCumulative -= m_currents[carrier_electrons][1];
}

void Cell::updateConcentrations()
{
	//electrons first, only in the film
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_concentrations[carrier_electrons][i] += m_currents[carrier_electrons][i] - m_currents[carrier_electrons][i + 1];
	}

	//cations 
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_cations][i] += m_currents[carrier_cations][i] - m_currents[carrier_cations][i + 1];
	}
	m_concentrations[carrier_cations][m_referencePoint] = m_saltConcentration;

	//anions
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_anions][i] += m_currents[carrier_anions][i] - m_currents[carrier_anions][i + 1];
	}
	m_concentrations[carrier_anions][m_referencePoint] = m_saltConcentration;

}

void Cell::loadState() //Kind of legacy function to load an existing state and continue the sim from there. Could probably be quite useful for analysis or to save time, but I never ended up using it due to worries about stability.
//e.g. it could be used to adjust a parameter mid-sim. (save state, change config, reload state, run).
{
	std::ifstream cFile("saveState.csv");
	if (cFile.is_open())
	{
		for (array_type::size_type i{ 0 }; i < m_size; ++i)
			cFile >> m_concentrations[carrier_electrons][i];

		for (array_type::size_type i{ 0 }; i < m_size; ++i)
			cFile >> m_concentrations[carrier_cations][i];

		for (array_type::size_type i{ 0 }; i < m_size; ++i)
			cFile >> m_concentrations[carrier_anions][i];
	}
	else
		throw - 1;
}

void Cell::resetInjection() //Just to make sure the injection function functions properly and does not run into negative concentration or something weird. The concentration are always updated immediately afterwards.
{
	if (m_concentrations[carrier_electrons][0] < 1)
		m_concentrations[carrier_electrons][0] = 1;
	if (m_concentrations[carrier_electrons][1] < 1)
		m_concentrations[carrier_electrons][1] = 1;
}

void Cell::changeBias(double vBiasChange) //legacy function, use the operators instead
{
	m_appliedBias += vBiasChange;
	m_electrostatic[es_potential][0] += vBiasChange;
} //increases the applied bias by vbiaschange

Cell& Cell::operator++()
{
	//manage the applied bias by incrementing (potentialArray is the big one, appliedBias is also incremented just to speed up convergense in calculatePotentialProfile)
	m_appliedBias += m_voltageIncrement;
	m_electrostatic[es_potential][0] += m_voltageIncrement;
	return *this;
}

Cell& Cell::operator--()
{
	//manage the applied bias by incrementing (potentialArray is the big one, appliedBias is also incremented just to speed up convergense in calculatePotentialProfile)
	m_appliedBias -= m_voltageIncrement;
	m_electrostatic[es_potential][0] -= m_voltageIncrement;
	return *this;
}

void Cell::midSave(std::ofstream& midf) //Save the current state of the sim as a .csv, newline splitted.
{
	{
		for (int i{ 0 }; i < 3; ++i)
		{
			for (array_type::size_type j{ 0 }; j < m_size; ++j)
			{
				midf << m_concentrations[i][j] << '\n';
			}
		}

		for (int i{ 0 }; i < 3; ++i)
		{
			for (array_type::size_type j{ 0 }; j < m_size; ++j)
			{
				midf << m_electrostatic[i][j] << '\n';
			}
		}

		for (int i{ 0 }; i < 3; ++i)
		{
			for (array_type::size_type j{ 0 }; j < m_size; ++j)
			{
				midf << m_currents[i][j] << '\n';
			}
		}
	}
}
