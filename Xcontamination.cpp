#include <cmath>
#include <cstddef>
#include <vector>
#include <array>
#include <numeric>
#include <iostream>
#include <fstream>
#include "Config.h"
#include "Electrochemistry.h" 
#include "Xcontamination.h"

using array_type = std::vector<std::vector<double>>;

Xcontamination::Xcontamination(settings_array& settings)
	:Cell{settings}, 
	m_currentConstantContaminant{ settings[s_contaminantMobility] * settings[s_dt] / settings[s_cellThickness] * m_size }
{
	m_contaminations = array_type(3, std::vector<double>(m_size));		//create the contaminant concentration array
	m_contaminationCurrents = array_type(3, std::vector<double>(m_size));	//create the contaminant current array
}

void Xcontamination::initializeConcentrations(double saltConcentration, double contaminantConcentration)
{
	//fills the concentration array with cations and anions
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_concentrations[carrier_cations][i] = saltConcentration * m_QDFillFactor;
		m_concentrations[carrier_anions][i] = saltConcentration * m_QDFillFactor;
	}

	for (array_type::size_type i{ m_interfacePoint }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_cations][i] = saltConcentration;
		m_concentrations[carrier_anions][i] = saltConcentration;
	}

	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_contaminations[contamination_X][i] = contaminantConcentration;
	}
}
void Xcontamination::calculatePotentialProfile()
{
	//This function takes in the reference to ePotentialArray, where the electrostatic potential profile and its 2 space derivates reside
	//It updates the array directly based on the input and thus does not need to return any values. vBias is also updated live to speed up convergence in the next steps

	//calculate the second space derivative of the electrostatic potential for every cell based on the Poisson equation
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantFilm * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i] 
																				- m_concentrations[carrier_anions][i] + m_contaminations[contamination_Xplus][i]);
	}

	for (array_type::size_type i{ m_interfacePoint }; i < (m_size - 1); ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantSolution * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i] 
																					- m_concentrations[carrier_anions][i] + m_contaminations[contamination_Xplus][i]);
	}

	double potAtReference{ 1.0 };
	while (std::abs(potAtReference) > 0.0001) //typically takes 5 rounds, quite a lot but for now its ok
	{
		//electricFieldStart is the initial guess for a boundary condition because the real boundary conditions are impossible to apply directly (assumes that the potential drops linearly over the cell at the start)
		//This boundary condition will then be updated based on the potential profile that is calculated
		//This process is repeated untal a potential profile is generated that satisfies the real boundary conditions (pot at negative electrode == applied potential && pot at reference electrode == 0)
		double electricFieldStart{ m_appliedBias / m_thickness };
		double differenceEF{ 1.0 };
		double electricFieldIntegral{};
		while (differenceEF > 0.00000001) //typically takes 2 round, very good!
		{
			m_electrostatic[es_electricField][0] = electricFieldStart; //apply guess bounday condition
			for (array_type::size_type i{ 0 }; i < (m_size - 2); ++i) //build electric field sequentially using said boundary condition
			{
				m_electrostatic[es_electricField][i + 1] = m_electrostatic[es_electricField][i] - m_electrostatic[es_secondDerivative][i + 1] * m_dx;
			}

			//check to see if real boundary condition is met (pot at negative elctrode == applied potential)
			electricFieldIntegral = std::accumulate(m_electrostatic[es_electricField].begin(), m_electrostatic[es_electricField].end(), 0.0) * m_dx;
			differenceEF = (m_appliedBias - electricFieldIntegral) / m_appliedBias;
			//and apply the updated boundary condition for next round
			electricFieldStart = electricFieldStart * (differenceEF + 1);
		}

		for (array_type::size_type i{ 0 }; i < (m_size - 1); ++i) // update potential profile using the aquired electrid field
		{
			m_electrostatic[es_potential][i + 1] = m_electrostatic[es_potential][i] - m_electrostatic[es_electricField][i] * m_dx;
		}

		//check to see if real boundary condition is met (pot at reference elctrode == 0)
		potAtReference = m_electrostatic[es_potential][m_referencePoint];
		//and update the bias to reflect this for the next round
		//divide by relatice position of the refernce point (default 1.5) for fastest conversion
		m_appliedBias += potAtReference / m_referencePositionRelative;
	}
}
inline double Xcontamination::negativeCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (-concentrationLeft * electricField + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
}

inline double Xcontamination::negativeCurrente(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	//if (concentrationLeft > 1e25)
	return (-concentrationLeft * (electricField - m_electronEnergyFactor * (concentrationLeft - concentrationRight) / pow((concentrationLeft > 0) ? concentrationLeft : 1, 0.3333333333)) + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
	//else
	//	return (-concentrationLeft * electricField  + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
}

inline double Xcontamination::positiveCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField)
{
	//returns the amount of positively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (concentrationRight * electricField + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
}

inline double Xcontamination::neutralCurrent(const double concentrationLeft, const double concentrationRight, const double curCon)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
}
void Xcontamination::react()
{

	//first we need to do the reaction of X to X- since it happens instantly when electrons are present
	for (array_type::size_type i{ 1 }; i < m_interfacePoint + 1; ++i)
	{
		//returns min(n,x) the lower one between electron and x concentrations in a certain cell
		m_contaminationCurrents[contamination_Xswitch][i] = (m_concentrations[carrier_electrons][i] - m_contaminations[contamination_Xplus][i]);
	}
/*
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		//update all the concentrations based on the oxidation reaction
		m_concentrations[carrier_electrons][i] -= m_contaminationCurrents[contamination_Xswitch][i];
	}
*/
	for (array_type::size_type i{ 1 }; i < m_interfacePoint+1; ++i)
	{
		//update all the concentrations based on the oxidation reaction
		m_concentrations[carrier_cations][i] -= m_contaminationCurrents[contamination_Xswitch][i];
	}

	for (array_type::size_type i{ 1 }; i < m_interfacePoint+1; ++i)
	{
		//update all the concentrations based on the oxidation reaction
		m_contaminations[contamination_Xplus][i] += m_contaminationCurrents[contamination_Xswitch][i];
	}
}
void Xcontamination::calculateCurrents()
{
	react(); //oxidation takes precedence over currents since it is a faster process

	//electrons first, only up to the interface, using the current function to fill the current array
	//the first cell is the electrode, see updateConcentrations on the exact working of the current there
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_currents[carrier_electrons][i] = negativeCurrent(m_concentrations[carrier_electrons][i - 1], m_concentrations[carrier_electrons][i],
			m_currentConstantElectrons, m_electrostatic[es_electricField][i - 1]);
	}

	//cations next, needs two different loops for the film and solution sections. Ions cannot enter the electrodes, so we dont need to consider the first and last cell
	for (array_type::size_type i{ 2 }; i < m_interfacePoint; ++i)
	{
		m_currents[carrier_cations][i] = positiveCurrent(m_concentrations[carrier_cations][i - 1], m_concentrations[carrier_cations][i],
			m_currentConstantCationsFilm, m_electrostatic[es_electricField][i - 1]);
	}

	for (array_type::size_type i{ m_interfacePoint + 1 }; i < (m_size - 1); ++i)
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
	m_currents[carrier_cations][m_interfacePoint] = positiveCurrent(m_concentrations[carrier_cations][m_interfacePoint - 1], m_concentrations[carrier_cations][m_interfacePoint] * m_QDFillFactor,
		m_currentConstantCationsFilm, m_electrostatic[es_electricField][m_interfacePoint - 1]);
	m_currents[carrier_anions][m_interfacePoint] = negativeCurrent(m_concentrations[carrier_anions][m_interfacePoint - 1], m_concentrations[carrier_anions][m_interfacePoint] * m_QDFillFactor,
		m_currentConstantAnionsFilm, m_electrostatic[es_electricField][m_interfacePoint - 1]);
/*
	//Added these for X and Xmin
	for (array_type::size_type i{ 2 }; i < (m_size - 1); ++i)
	{
		m_contaminationCurrents[contamination_X][i] = neutralCurrent(m_contaminations[contamination_X][i - 1], m_contaminations[contamination_X][i],
			m_currentConstantContaminant);
	}
*/
	for (array_type::size_type i{ 2 }; i < m_interfacePoint-1; ++i)
	{
		m_contaminationCurrents[contamination_Xplus][i] = positiveCurrent(m_contaminations[contamination_Xplus][i - 1], m_contaminations[contamination_Xplus][i],
			m_currentConstantContaminant, m_electrostatic[es_electricField][i - 1]);
	}

	m_currentCumulative -= m_currents[carrier_electrons][1];
}
void Xcontamination::updateConcentrations()
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
/*
	//X
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_contaminations[contamination_X][i] += m_contaminationCurrents[contamination_X][i] - m_contaminationCurrents[contamination_X][i + 1];
	}
*/
	//Xmin
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_contaminations[contamination_Xplus][i] += m_contaminationCurrents[contamination_Xplus][i] - m_contaminationCurrents[contamination_Xplus][i + 1];
	}
}
void Xcontamination::midSave(std::ofstream& midf)
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

		for (int i{ 0 }; i < 3; ++i)
		{
			for (array_type::size_type j{ 0 }; j < m_size; ++j)
			{
				midf << m_contaminations[i][j] << '\n';
			}
		}

		for (int i{ 0 }; i < 3; ++i)
		{
			for (array_type::size_type j{ 0 }; j < m_size; ++j)
			{
				midf << m_contaminationCurrents[i][j] << '\n';
			}
		}
	}
}

//from RI-reaction: contaminants using BV-kintetics
/*
RIreaction::RIreaction(settings_array& settings)
	:Cell{ settings },
	m_contaminantConcentration{ settings[s_sontaminantConcentration] },
	m_currentConstantX{ settings[s_contaminantMobility] * settings[s_dt] / settings[s_cellThickness] * m_size },
	m_k1forward{ settings[s_k1] * settings[s_dt] },
	m_k1backward{ settings[s_k1] * settings[s_dt] * exp(-settings[s_E1] * phys::q / phys::k / settings[s_temperature]) },
	m_k2{ settings[s_k2] * settings[s_dt] },
	m_maxTrapConcentration{ settings[s_maxTrapConcentration] }

{
	m_contaminants = array_type(3, std::vector<double>(m_size));		//create the special electrons concentration array
	m_contaminantCurrents = array_type(3, std::vector<double>(m_size));	//create the reation current array
}

void RIreaction::initializeConcentrations(double contaminantConcentration)
{
	//fills the concentration array with cations and anions
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_concentrations[carrier_cations][i] = m_saltConcentration * m_QDFillFactor;
		m_concentrations[carrier_anions][i] = m_saltConcentration * m_QDFillFactor;
		m_contaminants[contaminant_X][i] = contaminantConcentration;
	}

	for (array_type::size_type i{ m_interfacePoint }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_cations][i] = m_saltConcentration;
		m_concentrations[carrier_anions][i] = m_saltConcentration;
		m_contaminants[contaminant_X][i] = contaminantConcentration;
	}

	m_concentrations[carrier_electrons][0] = 1;
	m_concentrations[carrier_electrons][1] = 1;
}
void RIreaction::injectElectrons(const DOS_array& DOS)
{
	//These lines use the real DOS and Boltzmann-distribution to calculate the electron concentration.
		//DOS array also holds the "official" LUMO energy and energy integrating step dE
	double n{};
	double energy{ DOS[0] };
	double dE{ DOS[1] };
	//integrating from Estart to the end of DOS eV using steps of 0.01 eV
	for (DOS_array::size_type i{ 2 }; i < 300; ++i)
	{
		n += DOS[i] / (exp((energy - (m_negativeElectrodeWF - m_electrostatic[es_potential][0] + m_electrostatic[es_potential][1])) / m_energyConvert) + 1) * dE;
		energy += dE;
	}
	//m_concentrations[carrier_electrons][0] = n / 15;
	m_concentrations[carrier_electrons][0] = m_concentrations[carrier_electrons][0] * n / m_concentrations[carrier_electrons][1];
}
void RIreaction::calculatePotentialProfile()
{
	//This function takes in the reference to ePotentialArray, where the electrostatic potential profile and its 2 space derivates reside
	//It updates the array directly based on the input and thus does not need to return any values. vBias is also updated live to speed up convergence in the next steps

	//calculate the second space derivative of the electrostatic potential for every cell based on the Poisson equation
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantFilm * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i]
			- m_concentrations[carrier_anions][i] - m_contaminants[contaminant_Xmin][i]);
	}

	for (array_type::size_type i{ m_interfacePoint }; i < (m_size - 1); ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantSolution * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i]
			- m_concentrations[carrier_anions][i] - m_contaminants[contaminant_Xmin][i]);
	}

	double potAtReference{ 1.0 };
	while (std::abs(potAtReference) > 0.0001) //typically takes 2 rounds
	{
		//electricFieldStart is the initial guess for a boundary condition because the real boundary conditions are impossible to apply directly (assumes that the potential drops linearly over the cell at the start)
		//This boundary condition will then be updated based on the potential profile that is calculated
		//This process is repeated untal a potential profile is generated that satisfies the real boundary conditions (pot at negative electrode == applied potential && pot at reference electrode == 0)
		double electricFieldStart{ m_appliedBias / m_thickness };
		double differenceEF{ 1.0 };
		double electricFieldIntegral{};
		while (std::abs(differenceEF) > 0.00000001) //typically takes 2 round, very good!
		{
			m_electrostatic[es_electricField][0] = electricFieldStart; //apply guess bounday condition
			for (array_type::size_type i{ 0 }; i < (m_size - 2); ++i) //build electric field sequentially using said boundary condition
			{
				m_electrostatic[es_electricField][i + 1] = m_electrostatic[es_electricField][i] - m_electrostatic[es_secondDerivative][i + 1] * m_dx;
			}

			//check to see if real boundary condition is met (pot at negative elctrode == applied potential)
			electricFieldIntegral = std::accumulate(m_electrostatic[es_electricField].begin(), m_electrostatic[es_electricField].end(), 0.0) * m_dx;
			differenceEF = (m_appliedBias - electricFieldIntegral) / m_appliedBias;
			//and apply the updated boundary condition for next round
			electricFieldStart = electricFieldStart * (differenceEF + 1);
		}

		for (array_type::size_type i{ 0 }; i < (m_size - 1); ++i) // update potential profile using the aquired electrid field
		{
			m_electrostatic[es_potential][i + 1] = m_electrostatic[es_potential][i] - m_electrostatic[es_electricField][i] * m_dx;
		}

		//check to see if real boundary condition is met (pot at reference elctrode == 0)
		potAtReference = m_electrostatic[es_potential][m_referencePoint];
		//and update the bias to reflect this for the next round
		//divide by relatice position of the refernce point (default 1.5) for fastest conversion
		m_appliedBias += potAtReference / m_referencePositionRelative;
	}
}
void RIreaction::react()
{

	double kfi{ m_k1forward * exp(-(m_electrostatic[es_potential][0] - m_electrostatic[es_potential][1] - m_negativeElectrodeWF + m_LUMO) / m_energyConvert) };
	double kbi{ m_k1backward * m_k1forward / kfi };
	m_contaminantCurrents[switch_X][1] = (m_contaminants[contaminant_X][1] * (kfi < 1 ? kfi : 1)
		- m_contaminants[contaminant_Xmin][1] * (kbi < 1 ? kbi : 1));            //calculate the amount that reacts
	for (array_type::size_type i{ 2 }; i < m_interfacePoint + 1; ++i) //reaction only takes place in film
	{
		//
		double kf{ m_k1forward * exp(-(m_electrostatic[es_potential][i - 1] - m_electrostatic[es_potential][i]) / m_energyConvert) };
		double kb{ m_k1backward * m_k1forward / kf };
		//(m_maxTrapConcentration - m_electrons[electrons_trapped][i]) / m_maxTrapConcentration
		m_contaminantCurrents[switch_X][i] = (m_contaminants[contaminant_X][i] * (kf < 1 ? kf : 1)
			- m_contaminants[contaminant_Xmin][i] * (kb < 1 ? kb : 1));            //calculate the amount that reacts
		if (m_contaminantCurrents[switch_X][i] > m_concentrations[carrier_electrons][i - 1])
			m_contaminantCurrents[switch_X][i] = m_concentrations[carrier_electrons][i - 1];
	}
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_concentrations[carrier_electrons][i] -= m_contaminantCurrents[switch_X][i + 1];							//update the apporopriate concentrations
	}
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_contaminants[contaminant_X][i] -= m_contaminantCurrents[switch_X][i];							//update the apporopriate concentrations
	}
	for (array_type::size_type i{ 1 }; i < m_interfacePoint + 1; ++i)
	{
		m_contaminants[contaminant_Xmin][i] += m_contaminantCurrents[switch_X][i];
	}

	m_leakCurrentCumulative -= m_contaminantCurrents[contaminant_Xmin][m_size - 32];
	m_contaminants[contaminant_X][m_size - 32] += m_contaminants[contaminant_Xmin][m_size - 32];
	m_contaminants[contaminant_Xmin][m_size - 32] = 0;
}
inline double RIreaction::neutralCurrent(const double concentrationLeft, const double concentrationRight, const double curCon)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
}
void RIreaction::calculateCurrents()
{
	react(); //reaction takes precedence over currents (arbitrary but much easier this way)

	//electrons first, only up to the interface, using the current function to fill the current array
	//the first cell is the electrode, see updateConcentrations on the exact working of the current there
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_currents[carrier_electrons][i] = negativeCurrent(m_concentrations[carrier_electrons][i - 1], m_concentrations[carrier_electrons][i],
			m_currentConstantElectrons, m_electrostatic[es_electricField][i - 1]);
	}

	//cations next, needs two different loops for the film and solution sections. Ions cannot enter the electrodes, so we dont need to consider the first and last cell
	for (array_type::size_type i{ 2 }; i < m_interfacePoint; ++i)
	{
		m_currents[carrier_cations][i] = positiveCurrent(m_concentrations[carrier_cations][i - 1], m_concentrations[carrier_cations][i],
			m_currentConstantCationsFilm, m_electrostatic[es_electricField][i - 1]);
	}

	for (array_type::size_type i{ m_interfacePoint + 1 }; i < (m_size - 1); ++i)
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
	m_currents[carrier_cations][m_interfacePoint] = positiveCurrent(m_concentrations[carrier_cations][m_interfacePoint - 1], m_concentrations[carrier_cations][m_interfacePoint] * m_QDFillFactor,
		m_currentConstantCationsFilm, m_electrostatic[es_electricField][m_interfacePoint - 1]);
	m_currents[carrier_anions][m_interfacePoint] = negativeCurrent(m_concentrations[carrier_anions][m_interfacePoint - 1], m_concentrations[carrier_anions][m_interfacePoint] * m_QDFillFactor,
		m_currentConstantAnionsFilm, m_electrostatic[es_electricField][m_interfacePoint - 1]);

	//Added these for X and Xmin

	//for (array_type::size_type i{ 2 }; i < (m_size - 2); ++i)
	//{
	//	m_contaminantCurrents[contaminant_Xmin][i] = negativeCurrent(m_contaminants[contaminant_Xmin][i - 1], m_contaminants[contaminant_Xmin][i],
	//		m_currentConstantX, m_electrostatic[es_electricField][i - 1]);
	//}


	for (array_type::size_type i{ 2 }; i < (m_size - 1); ++i)
	{
		m_contaminantCurrents[contaminant_X][i] = neutralCurrent(m_contaminants[contaminant_X][i - 1], m_contaminants[contaminant_X][i],
			m_currentConstantX);
	}

	m_currentCumulative -= (m_currents[carrier_electrons][1] + m_contaminantCurrents[switch_X][1]);
}
void RIreaction::updateConcentrations()
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

	//Xmin
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_contaminants[contaminant_Xmin][i] += m_contaminantCurrents[contaminant_Xmin][i] - m_contaminantCurrents[contaminant_Xmin][i + 1];
	}

	//X
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_contaminants[contaminant_X][i] += m_contaminantCurrents[contaminant_X][i] - m_contaminantCurrents[contaminant_X][i + 1];
	}
	m_contaminants[contaminant_X][m_referencePoint] = m_contaminantConcentration;
}
void RIreaction::midSave(std::ofstream& midf)
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

		for (int i{ 0 }; i < 3; ++i)
		{
			for (array_type::size_type j{ 0 }; j < m_size; ++j)
			{
				midf << m_contaminants[i][j] << '\n';
			}
		}

		for (int i{ 0 }; i < 3; ++i)
		{
			for (array_type::size_type j{ 0 }; j < m_size; ++j)
			{
				midf << m_contaminantCurrents[i][j] << '\n';
			}
		}
	}
}
*///