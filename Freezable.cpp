#include <cmath>
#include <cstddef>
#include <vector>
#include <array>
#include <numeric>
#include <limits>
#include <iostream>
#include <fstream>
#include "Config.h"
#include "Electrochemistry.h" 
#include "Freezable.h"

double getDouble()
{
    while (true) // Loop until user enters a valid input
    {
        double x{};
        std::cin >> x;

        if (std::cin.fail()) // has a previous extraction failed?
        {
            // yep, so let's handle the failure
            std::cin.clear(); // put us back in 'normal' operation mode
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // and remove the bad input
        }
        else // else our extraction succeeded
        {
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            return x; // so return the value we extracted
        }
    }
}

Freezable::Freezable(Cell& cell)
	:Cell{ cell }
{
	m_frozencations = std::vector<double>(m_size);
	m_frozenanions = std::vector<double>(m_size);
}

void Freezable::calculatePotentialProfile()
{
	//This function takes in the reference to ePotentialArray, where the electrostatic potential profile and its 2 space derivates reside
	//It updates the array directly based on the input and thus does not need to return any values. vBias is also updated live to speed up convergence in the next steps

	//calculate the second space derivative of the electrostatic potential for every cell based on the Poisson equation
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantFilm * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i] + m_frozencations[i] 
																				- m_concentrations[carrier_anions][i] - m_frozenanions[i]);
	}

	for (array_type::size_type i{ m_interfacePoint }; i < (m_size - 1); ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantSolution * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i] + m_frozencations[i] - 
																				m_concentrations[carrier_anions][i] - m_frozenanions[i]);
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

void Freezable::updateConcentrations()
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

	//anions
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_anions][i] += m_currents[carrier_anions][i] - m_currents[carrier_anions][i + 1];
	}
}

void Freezable::setIonMobility(double dt)
{
    std::cout << "Enter the new ion mobility: ";
    double ionMobility{ getDouble() };
	m_currentConstantCationsFilm = ionMobility * dt / m_thickness * m_size;
	m_currentConstantAnionsFilm = ionMobility * dt / m_thickness * m_size;
}

void Freezable::fixIons() //fixes the fraction of ions that may not move again
{
    std::cout << "Enter the fraction of ions that should be frozen: ";
    double fraction{ getDouble() };
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		double totalcations{ m_concentrations[carrier_cations][i] + m_frozencations[i] };
		m_frozencations[i] = totalcations * fraction;
		m_concentrations[carrier_cations][i] = totalcations * (1-fraction);

		double totalanions{ m_concentrations[carrier_anions][i] + m_frozenanions[i] };
		m_frozenanions[i] = totalanions * fraction;
		m_concentrations[carrier_anions][i] = totalanions * (1 - fraction);
	}
}