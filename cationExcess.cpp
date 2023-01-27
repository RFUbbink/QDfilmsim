#include <cmath>
#include <cstddef>
#include <vector>
#include <array>
#include <numeric>
#include <iostream>
#include <fstream>
#include "Config.h"
#include "CationExcess.h"

using array_type = std::vector<std::vector<double>>;

CationExcess::CationExcess(settings_array& settings)
	:Cell{ settings }
	{
	m_cations = array_type(3, std::vector<double>(m_size));		//create the special electrons concentration array
	}

void CationExcess::calculatePotentialProfile()
{
	//This function takes in the reference to ePotentialArray, where the electrostatic potential profile and its 2 space derivates reside
	//It updates the array directly based on the input and thus does not need to return any values. vBias is also updated live to speed up convergence in the next steps

	//calculate the second space derivative of the electrostatic potential for every cell based on the Poisson equation
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantFilm * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i] + m_cations[cations_excess][i] - m_concentrations[carrier_anions][i]);
	}

	for (array_type::size_type i{ m_interfacePoint }; i < (m_size - 1); ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantSolution * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i] + m_cations[cations_excess][i] - m_concentrations[carrier_anions][i]);
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
		while (std::abs(differenceEF) > 0.00000001) //typically takes 2 rounds, very good!
		{
			m_electrostatic[es_electricField][0] = electricFieldStart; //apply guess boundary condition
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

		//check to see if real boundary condition is met (pot at reference electrode == 0)
		potAtReference = m_electrostatic[es_potential][m_referencePoint];
		//and update the bias to reflect this for the next round
		//divide by relatice position of the refernce point (default 1.5) for fastest conversion
		m_appliedBias += potAtReference / m_referencePositionRelative;
	}
}
void CationExcess::react()
{
	//first calculate how many cations are getting intercalated this step
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_cations[cations_toExcess][i] = m_concentrations[carrier_electrons][i] - m_cations[cations_excess][i];
	}

	//then update excess cation concetration
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_cations[cations_excess][i] += m_cations[cations_toExcess][i];
	}

	//then update cation and electron concetration
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_concentrations[carrier_cations][i] -= m_cations[cations_toExcess][i];
	}

	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_concentrations[carrier_electrons][i] -= m_cations[cations_toExcess][i];
	}

}
void CationExcess::calculateCurrents()
{
	react();
	//electrons first, only up to the interface, using the current function to fill the current array
	//the first step is the injection current

	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_currents[carrier_electrons][i] = negativeCurrent(m_concentrations[carrier_electrons][i - 1], m_concentrations[carrier_electrons][i],
			m_currentConstantElectrons, m_electrostatic[es_electricField][i - 1]);
	}


	//cations next, needs two different loops for the film and solution sections. Ions cannot enter the electrodes, so we dont need to consider the first and last cell
	//also, we use the maxed out current in the film because this is the only place it can be relevant
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

	//Excess cations
	for (array_type::size_type i{ 2 }; i < m_interfacePoint; ++i)
	{
		m_cations[cations_current][i] = positiveCurrent(m_cations[cations_excess][i - 1], m_cations[cations_excess][i],
			m_currentConstantCationsFilm, m_electrostatic[es_electricField][i - 1]);
	}

	for (array_type::size_type i{ m_interfacePoint + 1 }; i < (m_size - 1); ++i)
	{
		m_cations[cations_current][i] = positiveCurrent(m_cations[cations_excess][i - 1], m_cations[cations_excess][i],
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
	m_cations[cations_current][m_interfacePoint] = positiveCurrent(m_cations[cations_excess][m_interfacePoint - 1], m_cations[cations_excess][m_interfacePoint] * m_QDFillFactor,
		m_currentConstantCationsFilm, m_electrostatic[es_electricField][m_interfacePoint - 1]);

	//update total electrons that have entered since the last getCurrent() call:
	m_currentCumulative -= m_currents[carrier_electrons][1] * m_currentConvert;

}
void CationExcess::updateConcentrations()
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

	//excess cations 
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_cations[cations_excess][i] += m_cations[cations_current][i] - m_cations[cations_current][i + 1];
	}

	//anions
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_anions][i] += m_currents[carrier_anions][i] - m_currents[carrier_anions][i + 1];
	}
	m_concentrations[carrier_anions][m_referencePoint] = m_saltConcentration;

}
void CationExcess::midSave(std::ofstream& midf)
{
	{
		//concentrationArray
		for (int i{ 0 }; i < 3; ++i)
		{
			for (array_type::size_type j{ 0 }; j < m_size; ++j)
			{
				midf << m_concentrations[i][j] << '\n';
			}
		}

		//electrostaticArray
		for (int i{ 0 }; i < 3; ++i)
		{
			for (array_type::size_type j{ 0 }; j < m_size; ++j)
			{
				midf << m_electrostatic[i][j] << '\n';
			}
		}

		//currentArray
		for (int i{ 0 }; i < 3; ++i)
		{
			for (array_type::size_type j{ 0 }; j < m_size; ++j)
			{
				midf << m_currents[i][j] << '\n';
			}
		}

		//interlatated cations
		for (array_type::size_type j{ 0 }; j < m_size; ++j)
		{
			midf << m_cations[cations_excess][j] << '\n';
		}
	}
}