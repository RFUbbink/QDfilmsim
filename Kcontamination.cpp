#include <cmath>
#include <cstddef>
#include <vector>
#include <array>
#include <numeric>
#include <iostream>
#include <fstream>
#include "Config.h"
#include "Electrochemistry.h" 
#include "Kcontamination.h"

using array_type = std::vector<std::vector<double>>;

Kcontamination::Kcontamination(settings_array& settings)
	:Cell{ settings },
	m_currentConstantContaminant{ settings[s_contaminantMobility] * settings[s_dt] / settings[s_cellThickness] * m_size },
	m_maximumElectronConcentration{ settings[s_maximumElectronConcentration] }
{
	m_contaminations = array_type(3, std::vector<double>(m_size));		//create the contaminant concetration array
	m_contaminationCurrents = array_type(3, std::vector<double>(m_size));	//create the contaminant current array
}

void Kcontamination::injectElectrons()
{
	double inject{ (m_densityOfStates * exp((m_electrostatic[es_potential][1] - m_electrostatic[es_potential][0] - m_injectionBarrier) / m_energyConvert)) };
	m_concentrations[0][0] = ((inject < m_maximumElectronConcentration) ? inject : m_maximumElectronConcentration);
}
void Kcontamination::initializeConcentrations(double saltConcentration, double contaminantConcentration)
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
		m_contaminations[contamination_K][i] = contaminantConcentration;
	}

	for (array_type::size_type i{ 0 }; i < m_size; ++i)
	{
		m_contaminations[contamination_switch_maxCon][i] = m_maximumElectronConcentration;
	}
}
void Kcontamination::calculatePotentialProfile()
{
	//This function takes in the reference to ePotentialArray, where the electrostatic potential profile and its 2 space derivates reside
	//It updates the array directly based on the input and thus does not need to return any values. vBias is also updated live to speed up convergence in the next steps

	//calculate the second space derivative of the electrostatic potential for every cell based on the Poisson equation
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantFilm * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i]
			- m_concentrations[carrier_anions][i] - m_contaminations[contamination_QDKmin][i]);
	}

	for (array_type::size_type i{ m_interfacePoint }; i < (m_size - 1); ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantSolution * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i]
			- m_concentrations[carrier_anions][i]);
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
inline double Kcontamination::negativeCurrentmax(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double maximum)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	double current{ (-concentrationLeft * electricField + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon };
	return ((current < maximum) ? current : maximum);
}
inline double Kcontamination::neutralCurrent(const double concentrationLeft, const double concentrationRight, const double curCon)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
}
void Kcontamination::react()
{
	//first we need to do the reaction of K to QDK- since it happens instantly when electrons are present
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		//returns min(n,x) the lower one between electron and K concentrations in a certain cell
		m_contaminationCurrents[contamination_switch_maxCon][i] = ((m_concentrations[carrier_electrons][i] < m_contaminations[contamination_K][i])
			? m_concentrations[carrier_electrons][i] : m_contaminations[contamination_K][i]);
	}

	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		//update all the concentrations based on the oxidation reaction
		m_concentrations[carrier_electrons][i] -= m_contaminationCurrents[contamination_switch_maxCon][i];
	}

	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		//update all the concentrations based on the oxidation reaction
		m_contaminations[contamination_K][i] -= m_contaminationCurrents[contamination_switch_maxCon][i];
	}

	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		//update all the concentrations based on the oxidation reaction
		m_contaminations[contamination_QDKmin][i] += m_contaminationCurrents[contamination_switch_maxCon][i];
	}

	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		//update all the concentrations based on the oxidation reaction
		m_contaminations[contamination_switch_maxCon][i] -= m_contaminationCurrents[contamination_switch_maxCon][i];
	}

}
void Kcontamination::calculateCurrents()
{
	react(); //oxidation takes precedence over currents since it is a faster process

	//electrons first, only up to the interface, using the current function to fill the current array
	//the first cell is the electrode, see updateConcentrations on the exact working of the current there
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_currents[carrier_electrons][i] = negativeCurrentmax(m_concentrations[carrier_electrons][i - 1], m_concentrations[carrier_electrons][i],
			m_currentConstantElectrons, m_electrostatic[es_electricField][i - 1], (m_contaminations[contamination_switch_maxCon][i]-m_concentrations[carrier_electrons][i]));
	}

	//cations next, needs two different loops for the film and solution sections. Ions cannot enter the electrodes, so we dont need to consider the first and last cell
	for (array_type::size_type i{ 2 }; i < m_interfacePoint; ++i)
	{
		m_currents[carrier_cations][i] = positiveCurrentmax(m_concentrations[carrier_cations][i - 1], m_concentrations[carrier_cations][i],
			m_currentConstantCationsFilm, m_electrostatic[es_electricField][i - 1], (m_maximumIonConcentration - m_concentrations[carrier_cations][i - 1]));
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
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_contaminationCurrents[contamination_K][i] = neutralCurrent(m_contaminations[contamination_K][i - 1], m_contaminations[contamination_K][i],
			m_currentConstantContaminant);
	}


}
void Kcontamination::updateConcentrations()
{
	//electrons first, only in the film
	for (array_type::size_type i{ 0 }; i < m_interfacePoint; ++i)
	{
		m_concentrations[carrier_electrons][i] += m_currents[carrier_electrons][i] - m_currents[carrier_electrons][i + 1];
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

	//K
	for (array_type::size_type i{ 2 }; i < (m_size - 2); ++i)
	{
		m_contaminations[contamination_K][i] += m_contaminationCurrents[contamination_K][i] - m_contaminationCurrents[contamination_K][i + 1];
	}


}
void Kcontamination::midSave(std::ofstream& midf)
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