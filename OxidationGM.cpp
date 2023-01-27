#include <cmath>
#include <cstddef>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include "Config.h"
#include "Electrochemistry.h" 
#include "OxidationGM.h"

using array_type = std::vector<std::vector<double>>;

OxidationGM::OxidationGM(settings_array& settings, DOS_array& RR)
	:Cell{ settings },
	m_contaminantConcentration{ settings[s_sontaminantConcentration] },
	m_currentConstantX{ settings[s_contaminantMobility] * settings[s_dt] / settings[s_cellThickness] * m_size },
	m_k{settings[s_k1] * settings[s_dt]},
	m_RRlookup{RR},
	m_stepCounter{ 0 }

{
	m_contaminants = array_type(3, std::vector<double>(m_size));		//create the oxidadant/reductant concentration array
	m_contaminantCurrents = array_type(3, std::vector<double>(m_size));	//create the reation current array
	m_reactionRates = array_type(3, std::vector<double>(m_size));		//create the array that holds the reaction rates
}

void OxidationGM::initializeConcentrations(double contaminantConcentration)
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

void OxidationGM::calculatePotentialProfile()
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

void OxidationGM::updateRates()
{
	DOS_array::size_type pos{};
	double n{};
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		n = m_concentrations[carrier_electrons][i];
		if (n > m_RRlookup[0])
		{
			std::cerr << "Electron concentration exceeded maximum! Adjust the reaction rate table to accomodate higher concentrations.";
			std::exit(1);
		}
		pos =  static_cast<DOS_array::size_type>(1 + n/m_RRlookup[0]*200);
		//this is the left position, pos + 1 is the first index where n is higher than the value we are searching for
		
		m_reactionRates[rate_forward][i] = m_k * (m_RRlookup[pos]+ (200*n/m_RRlookup[0]-(pos-1))*(m_RRlookup[pos+1] - m_RRlookup[pos]));
		m_reactionRates[rate_backward][i] = m_k * (m_RRlookup[pos + 200] + (200 * n / m_RRlookup[0] - (pos - 1))*(m_RRlookup[pos + 201] - m_RRlookup[pos + 200]));
	}
}

void OxidationGM::react()
{
	if (++m_stepCounter == 10) //Rates are updated every 10 steps. Can be adjusted here easily to every 5 or every step. 
	{
		updateRates();
		m_stepCounter = 0;
	}


	for (array_type::size_type i{ 2 }; i < m_interfacePoint + 1; ++i)
	{
		m_contaminantCurrents[switch_X][i] = m_contaminants[contaminant_X][i] * m_reactionRates[contaminant_X][i]
												- m_contaminants[contaminant_Xmin][i] * m_reactionRates[contaminant_Xmin][i]; //swtich_x is the amount that goes from X to Xmin
	}

	for (array_type::size_type i{ 2 }; i < m_interfacePoint + 1; ++i)
	{
		m_concentrations[carrier_electrons][i] -= m_contaminantCurrents[switch_X][i];
	}

	for (array_type::size_type i{ 2 }; i < m_interfacePoint + 1; ++i)
	{
		m_contaminants[contaminant_X][i] -= m_contaminantCurrents[switch_X][i];
	}

	for (array_type::size_type i{ 2 }; i < m_interfacePoint + 1; ++i)
	{
		m_contaminants[contaminant_Xmin][i] += m_contaminantCurrents[switch_X][i];
	}

	m_contaminants[contaminant_X][m_referencePoint - 1] += m_contaminants[contaminant_Xmin][m_referencePoint - 1];
	m_contaminants[contaminant_Xmin][m_referencePoint - 1] = 0;
}

inline double OxidationGM::negativeCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (-concentrationLeft * electricField + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
}

inline double OxidationGM::negativeCurrente(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	//if (concentrationLeft > 1e25)
	return (-concentrationLeft * (electricField - m_electronEnergyFactor * (concentrationLeft - concentrationRight) / pow((concentrationLeft > 0) ? concentrationLeft : 1, 0.3333333333)) + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
	//else
	//	return (-concentrationLeft * electricField  + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
}

inline double OxidationGM::positiveCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField)
{
	//returns the amount of positively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (concentrationRight * electricField + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
}

inline double OxidationGM::neutralCurrent(const double concentrationLeft, const double concentrationRight, const double curCon)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
}

void OxidationGM::calculateCurrents()
{
	react(); //reaction takes precedence over currents (arbitrary but much easier this way)

	//electrons first, only up to the interface, using the current function to fill the current array
	//the first cell is the electrode, see updateConcentrations on the exact working of the current there
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


	//for (array_type::size_type i{ 2 }; i < (m_size - 1); ++i)
	//{
	//	m_contaminantCurrents[contaminant_X][i] = neutralCurrent(m_contaminants[contaminant_X][i - 1], m_contaminants[contaminant_X][i],
	//		m_currentConstantX);
	//}

	m_currentCumulative -= (m_currents[carrier_electrons][1] + m_contaminantCurrents[switch_X][1]);
}

void OxidationGM::updateConcentrations()
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
	//for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	//{
	//	m_contaminants[contaminant_Xmin][i] += m_contaminantCurrents[contaminant_Xmin][i] - m_contaminantCurrents[contaminant_Xmin][i + 1];
	//}

	//X
	//for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	//{
	//	m_contaminants[contaminant_X][i] += m_contaminantCurrents[contaminant_X][i] - m_contaminantCurrents[contaminant_X][i + 1];
	//}
	//m_contaminants[contaminant_X][m_referencePoint] = m_contaminantConcentration;
}

void OxidationGM::midSave(std::ofstream& midf)
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
