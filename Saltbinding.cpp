#include <cmath>
#include <cstddef>
#include <vector>
#include <array>
#include <numeric>
#include <iostream>
#include <fstream>
#include "Config.h"
#include "Electrochemistry.h" 
#include "Saltbinding.h"

Saltbinding::Saltbinding(settings_array& settings)
	:Cell{ settings },
	m_captureCoefficient{1e-25},
	m_dissociationCoeffecient{3.02e-28},
	m_dt{settings[s_dt]}
{
	m_salt = array_type(3, std::vector<double>(m_size));		//create the potential array
}

void Saltbinding::react()
{
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_salt[salt_toSalt][i] = (m_captureCoefficient * m_concentrations[carrier_anions][i] * (m_concentrations[carrier_cations][i] - m_concentrations[carrier_electrons][i])
								- 1e27 * m_dissociationCoeffecient * m_salt[salt_salt][i])*m_dt;
	}

	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_salt[salt_salt][i] += m_salt[salt_toSalt][i];
	}

	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_cations][i] -= m_salt[salt_toSalt][i];
	}

	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_anions][i] -= m_salt[salt_toSalt][i];
	}
}
void Saltbinding::initializeConcentrations(double contaminantConcentration)
{
	contaminantConcentration;
	//fills the concentration array with cations and anions
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_concentrations[carrier_cations][i] = m_saltConcentration * m_QDFillFactor;
		m_concentrations[carrier_anions][i] = m_saltConcentration * m_QDFillFactor;
	}

	for (array_type::size_type i{ m_interfacePoint }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_cations][i] = m_saltConcentration;
		m_concentrations[carrier_anions][i] = m_saltConcentration;
	}
	
	do{react();} 
	while (m_salt[salt_toSalt][60] > m_salt[salt_salt][60] / 1000);

	m_saltConcentration = m_salt[salt_salt][60];

	m_concentrations[carrier_electrons][0] = 1;
	m_concentrations[carrier_electrons][1] = 1;
}
void Saltbinding::calculateCurrents()
{
	react();

	//electrons first, only up to the interface, using the current function to fill the current array
	//the first step is the injection current
/*
	m_currents[carrier_electrons][1] = negativeCurrent(m_concentrations[carrier_electrons][0], m_concentrations[carrier_electrons][1],
				m_currentConstantElectrons, m_electrostatic[es_electricField][0]);

	for (array_type::size_type i{ m_interfacePoint - 20 }; i < m_interfacePoint; ++i)
	{
		m_currents[carrier_electrons][i] = negativeCurrentmax(m_concentrations[carrier_electrons][i - 1], m_concentrations[carrier_electrons][i],
			m_currentConstantElectrons, m_electrostatic[es_electricField][i - 1], (1.5*m_concentrations[carrier_electrons][1]- m_concentrations[carrier_electrons][i]));
	}
*/

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


	//update total electrons that have entered since the last getCurrent() call:
	m_currentCumulative -= m_currents[carrier_electrons][1];

}
void Saltbinding::midSave(std::ofstream& midf)
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

	for (array_type::size_type j{ 0 }; j < m_size; ++j)
	{
		midf << m_salt[salt_salt][j] << '\n';
	}

}
