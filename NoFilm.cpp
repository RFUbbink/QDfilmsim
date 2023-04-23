#include <cmath>
#include <cstddef>
#include <vector>
#include <array>
#include <numeric>
#include <iostream>
#include <fstream>
#include "Config.h"
#include "NoFilm.h" 

using array_type = std::vector<std::vector<double>>;

//This class can be used to run a "normal" flat-electrode redox simulation (although ridiculously slow compared to other methods).
//You would simulate an electrochemical cell with ions present as electrolyte, and 1 redox species at a certain concentration.
//All species are in oxidized state at the start of simulation.
//Again, many functions are the same as in Electrochemistry, so look for comments there if unclear

NoFilmCell::NoFilmCell(settings_array& settings) //initialize all the constants necessary for operation, using settings array
//vBias is the TOTAL potential drop over the entire system == the voltage difference between the working and counter electrode
//an estimation of the initial voltage drop over the counterelectrode is made to calculate the initial Vb.
	:m_appliedBias{ settings[s_startVoltage] },
	m_voltageIncrement{ settings[s_voltageIncrement] },
	m_saltConcentration{ settings[s_ionConcentration] },
	m_Xconcentration {settings[s_redoxSpeciesConcentration]},
	m_size{ static_cast<array_type::size_type>(settings[s_amountOfCells]) },
	m_referencePoint{ static_cast<array_type::size_type>(m_size / 2) },
	m_referencePositionRelative{static_cast<double>(m_referencePoint -1) / static_cast<double>(m_size) },
	m_thickness{ settings[s_cellThickness] },
	m_dx{ settings[s_cellThickness] / m_size },
	m_E0{ settings[s_LUMO]},
	m_currentConstantX{ settings[s_electronMobility] * settings[s_dt] / m_dx },
	m_currentConstantCations{ settings[s_cationMobilitySolution] * settings[s_dt] / m_dx },
	m_currentConstantAnions{ settings[s_anionMobilitySolution] * settings[s_dt] / m_dx },
	m_energyConvert{ physics::k * settings[s_temperature] / physics::q },
	m_energyConvertx{ physics::k * settings[s_temperature] / physics::q / m_dx },
	m_poissonConstantSolution{ -physics::q / physics::eps0 / settings[s_epsilonrSolution] },
	m_currentConvert{ physics::q * m_dx / settings[s_dt] / 10000 },
	m_Nernst{ physics::F / physics::R / settings[s_temperature]}
{
	m_electrostatic = array_type(3, std::vector<double>(m_size));		//create the potential array
	m_electrostatic[es_potential][0] = settings[s_startVoltage];		//apply the inital bias: updated everytime there is a voltage increment.
																		//This array value should be constant unless there is an increment
	m_concentrations = array_type(4, std::vector<double>(m_size));		//create the concentration array
	m_currents = array_type(4, std::vector<double>(m_size));			//create the current array
}

double NoFilmCell::getCurrent() //returns the electron current at a given time (measured in the first cell next to the negative electrode)
{
	double curcur{ m_currentCumulative * m_currentConvert };
	m_currentCumulative = 0;
	return curcur;
}

double NoFilmCell::getLeakCurrent() //returns the electron current at a given time (measured in the first cell next to the negative electrode)
{
	return 0;
}

double NoFilmCell::getVoltage()
{
	return m_electrostatic[es_potential][0];
}

void NoFilmCell::resetInjection() {
	
}

void NoFilmCell::injectElectrons(const DOS_array& DOS)
{
	//This function actually performs the electrochemical reaction at the electrode.
	//Assumed is Nernstian equilibrium, reaction kinetics are neglected.
	//It is a modified injectElectrons() from the other classes, so they can be used as same Types in the main() and runIV() functions. 
	DOS;
	double totalX{ m_concentrations[carrier_X][1] + m_concentrations[carrier_Xmin][1] };
	double equilibriumRatio{ exp((m_E0 - m_electrostatic[es_potential][0] + m_electrostatic[es_potential][1]) * m_Nernst) };
	m_concentrations[carrier_Xmin][1] = totalX * equilibriumRatio / (1 + equilibriumRatio);
	m_currentCumulative -= m_concentrations[carrier_Xmin][1] - (totalX- m_concentrations[carrier_X][1]);
	m_concentrations[carrier_X][1] = totalX / (1 + equilibriumRatio); 
}

void NoFilmCell::calculatePotentialProfile()
{
	//This function takes in the reference to ePotentialArray, where the electrostatic potential profile and its 2 space derivates reside
	//It updates the array directly based on the input and thus does not need to return any values. vBias is also updated live to speed up convergence in the next steps

	//calculate the second space derivative of the electrostatic potential for every cell based on the Poisson equation

	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantSolution * (-m_concentrations[carrier_Xmin][i] + m_concentrations[carrier_cations][i] - m_concentrations[carrier_anions][i]);
	}
	//std::cout << "\nNewPPcall\n";
	double potAtReference{ 1.0 };
	while (std::abs(potAtReference) > 0.0001) //typically takes 5 rounds, quite a lot but for now its ok
	{
		//electricFieldStart is the initial guess for a boundary condition because the real boundary conditions are impossible to apply directly (assumes that the potential drops linearly over the cell at the start)
		//This boundary condition will then be updated based on the potential profile that is calculated
		//This process is repeated untal a potential profile is generated that satisfies the real boundary conditions (pot at negative electrode == applied potential && pot at reference electrode == 0)
		double electricFieldStart{ m_appliedBias / m_thickness * m_size/(m_size-1) };
		double differenceEF{ 1.0 };
		double electricFieldIntegral{};
		while (std::abs(differenceEF) > 0.00001) //typically takes 2 rounds, very good!
		{
			m_electrostatic[es_electricField][0] = electricFieldStart; //apply guess boundary condition
			for (array_type::size_type i{ 0 }; i < (m_size - 2); ++i) //build electric field sequentially using said boundary condition
			{
				m_electrostatic[es_electricField][i + 1] = m_electrostatic[es_electricField][i] - m_electrostatic[es_secondDerivative][i + 1] * m_dx;
			}

			//check to see if real boundary condition is met (pot at negative elctrode == applied potential)
			electricFieldIntegral = std::accumulate(m_electrostatic[es_electricField].begin(), m_electrostatic[es_electricField].end(), 0.0) * m_dx;
			differenceEF = (m_appliedBias - electricFieldIntegral) / m_appliedBias;
			//std::cout << differenceEF << '\n';
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
		//divide by relative position of the refernce point (default 1.5) for fastest conversion
		m_appliedBias += potAtReference / m_referencePositionRelative;
	}
}

void NoFilmCell::initializeConcentrations()
{

	//fills the concentration array with cations and anions
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_cations][i] = m_saltConcentration;
		m_concentrations[carrier_anions][i] = m_saltConcentration;
	}
	//Sets the concentration of redos species throughout the cell
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_X][i] = m_Xconcentration;
	}
}

inline double NoFilmCell::negativeCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (-concentrationLeft * electricField + eCon * (concentrationLeft - concentrationRight)) * curCon;
}

inline double NoFilmCell::neutralCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double eCon)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (eCon * (concentrationLeft - concentrationRight)) * curCon;
}

inline double NoFilmCell::positiveCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon)
{
	//returns the amount of positively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (concentrationRight * electricField + eCon * (concentrationLeft - concentrationRight)) * curCon;
}

void NoFilmCell::calculateCurrents()
{

	//Added these for X and Xmin
	//Xmin
	for (array_type::size_type i{ 2 }; i < (m_size - 1); ++i)
	{
		m_currents[carrier_Xmin][i] = negativeCurrent(m_concentrations[carrier_Xmin][i - 1], m_concentrations[carrier_Xmin][i],
			m_currentConstantX, m_electrostatic[es_electricField][i - 1], m_energyConvertx);
	}

	//X
	for (array_type::size_type i{ 2 }; i < (m_size - 1); ++i)
	{
		m_currents[carrier_X][i] = neutralCurrent(m_concentrations[carrier_X][i - 1], m_concentrations[carrier_X][i],
			m_currentConstantX, m_energyConvertx);
	}

	//cations next, needs two different loops for the film and solution sections. Ions cannot enter the electrodes, so we dont need to consider the first and last cell
	for (array_type::size_type i{2}; i < (m_size - 1); ++i)
	{
		m_currents[carrier_cations][i] = positiveCurrent(m_concentrations[carrier_cations][i - 1], m_concentrations[carrier_cations][i],
			m_currentConstantCations, m_electrostatic[es_electricField][i - 1], m_energyConvertx);
	}

	//the anions which are very similar to cations
	for (array_type::size_type i{ 2}; i < (m_size - 1); ++i)
	{
		m_currents[carrier_anions][i] = negativeCurrent(m_concentrations[carrier_anions][i - 1], m_concentrations[carrier_anions][i],
			m_currentConstantAnions, m_electrostatic[es_electricField][i - 1], m_energyConvertx);
	}
}

void NoFilmCell::updateConcentrations()
{
	
	//Xmin
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_Xmin][i] += m_currents[carrier_Xmin][i] - m_currents[carrier_Xmin][i + 1];
	}

	//X
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_X][i] += m_currents[carrier_X][i] - m_currents[carrier_X][i + 1];
	}
	//cations 
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_cations][i] += m_currents[carrier_cations][i] - m_currents[carrier_cations][i + 1];
	}

	//anions
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_anions][i] += m_currents[carrier_anions][i] - m_currents[carrier_anions][i + 1];
	}
}

NoFilmCell& NoFilmCell::operator++()
{
	//manage the applied bias by incrementing (potentialArray is the big one, appliedBias is also incrementen just to speed up convergense in calculatePotentialProfile)
	m_appliedBias += m_voltageIncrement;
	m_electrostatic[es_potential][0] += m_voltageIncrement;
	return *this;
}

NoFilmCell& NoFilmCell::operator--()
{
	//manage the applied bias by incrementing (potentialArray is the big one, appliedBias is also incrementen just to speed up convergense in calculatePotentialProfile)
	m_appliedBias -= m_voltageIncrement;
	m_electrostatic[es_potential][0] -= m_voltageIncrement;
	return *this;
}

void NoFilmCell::midSave(std::ofstream& midf)
{
	{
		for (int i{ 0 }; i < 4; ++i)
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

		for (int i{ 0 }; i < 4; ++i)
		{
			for (array_type::size_type j{ 0 }; j < m_size; ++j)
			{
				midf << m_currents[i][j] << '\n';
			}
		}
	}
}

