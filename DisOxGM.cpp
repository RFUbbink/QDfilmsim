#include <cmath>
#include <cstddef>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include "Config.h"
#include "Discontinuous.h" 
#include "DisOxGM.h"

using array_type = std::vector<std::vector<double>>;

DisOxGM::DisOxGM(settings_array& settings, DOS_array& RR)
	:DisCell{ settings },
	m_contaminantConcentration{ settings[s_sontaminantConcentration] },
	m_currentConstantXFilm{ settings[s_contaminantMobility] * settings[s_dt] / m_dxf },
	m_currentConstantXSolution1{ settings[s_contaminantMobility] * settings[s_dt] / m_dxs1 },
	m_currentConstantXSolution2{ settings[s_contaminantMobility] * settings[s_dt] / m_dxs2 },
	m_k{settings[s_k1]*settings[s_dt]},
	m_E0{settings[s_E1]},
	m_RRlookup{ RR },
	m_stepCounter{ 0 }

{
	m_contaminants = array_type(3, std::vector<double>(m_size));		//create the oxidadant/reductant concentration array
	m_contaminantCurrents = array_type(3, std::vector<double>(m_size));	//create the reation current array
	m_reactionRates = array_type(3, std::vector<double>(m_size));		//create the array that holds the reaction rates
	std::cout << "Currently applying Gerisher-Marcus kinetics\n";  //Switch m_k to settings[s_dt]
	//std::cout << "Currently applying Butler-Volmer kinetics (or trying anyway)\n"; //Switch m_k to settings[s_k1]*settings[s_dt]
	//std::cout << "Currently applying an 'Nc value' of 1e27 for electron rate\n";
	//std::cout << "Diffusion of conctaminants disabled!\n";
}

inline double DisOxGM::negativeCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (-concentrationLeft * electricField + eCon * (concentrationLeft - concentrationRight)) * curCon;
}

inline double DisOxGM::negativeCurrente(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	//if (concentrationLeft > 1e25)
	return (-concentrationLeft * (electricField - m_electronEnergyFactor * (concentrationLeft - concentrationRight) / pow((concentrationLeft > 0) ? concentrationLeft : 1, 0.3333333333)) + eCon * (concentrationLeft - concentrationRight)) * curCon;
	//else
	//	return (-concentrationLeft * electricField  + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
}

inline double DisOxGM::positiveCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon)
{
	//returns the amount of positively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (concentrationRight * electricField + eCon * (concentrationLeft - concentrationRight)) * curCon;
}

inline double DisOxGM::neutralCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double eCon)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (eCon * (concentrationLeft - concentrationRight)) * curCon;
}

void DisOxGM::initializeConcentrations(double contaminantConcentration)
{
	std::cout << "Electron averaging on!\nSlow mobility near interface on!\n";
	//fills the concentration array with cations and anions and the FILM ONLY with contaminant
	for (array_type::size_type i{ 1 }; i < m_interfacePoint + 1; ++i)
	{
		m_concentrations[carrier_cations][i] = m_saltConcentration * m_QDFillFactor;
		m_concentrations[carrier_anions][i] = m_saltConcentration * m_QDFillFactor;
		m_contaminants[contaminant_X][i] = contaminantConcentration;
	}

	for (array_type::size_type i{ m_interfacePoint + 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_cations][i] = m_saltConcentration;
		m_concentrations[carrier_anions][i] = m_saltConcentration;
		m_contaminants[contaminant_X][i] = contaminantConcentration;
	}

	m_concentrations[carrier_electrons][0] = 1;
	m_concentrations[carrier_electrons][1] = 1;
}

void DisOxGM::calculatePotentialProfile()
{
	//std::cout << m_dxf << '\n' << m_dxs1 << '\n' << m_dxs2 << '\n';
	//This function takes in the reference to ePotentialArray, where the electrostatic potential profile and its 2 space derivates reside
	//It updates the array directly based on the input and thus does not need to return any values. vBias is also updated live to speed up convergence in the next steps

	//calculate the second space derivative of the electrostatic potential for every cell based on the Poisson equation
	for (array_type::size_type i{ 1 }; i < m_interfacePoint + 1; ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantFilm * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i] 
																				- m_concentrations[carrier_anions][i] - m_contaminants[contaminant_Xmin][i]);
	}

	for (array_type::size_type i{ m_interfacePoint + 1 }; i < (m_size - 1); ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantSolution * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i] 
																				- m_concentrations[carrier_anions][i] - m_contaminants[contaminant_Xmin][i]);
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
		while (std::abs(differenceEF) > 0.00001) //typically takes 2 rounds, very good!
		{
			m_electrostatic[es_electricField][0] = electricFieldStart; //apply guess boundary condition
			for (array_type::size_type i{ 0 }; i < m_interfacePoint - 5; ++i) //build electric field sequentially using said boundary condition
			{
				m_electrostatic[es_electricField][i + 1] = m_electrostatic[es_electricField][i] - m_electrostatic[es_secondDerivative][i + 1] * m_dxf;
			}


			for (array_type::size_type i{ m_interfacePoint - 5 }; i < (m_interfacePoint + m_Isize); ++i) //build electric field sequentially using said boundary condition
			{
				m_electrostatic[es_electricField][i + 1] = m_electrostatic[es_electricField][i] - m_electrostatic[es_secondDerivative][i + 1] * m_dxs1;
			}

			for (array_type::size_type i{ m_interfacePoint + m_Isize }; i < (m_size - 2); ++i) //build electric field sequentially using said boundary condition
			{
				m_electrostatic[es_electricField][i + 1] = m_electrostatic[es_electricField][i] - m_electrostatic[es_secondDerivative][i + 1] * m_dxs2;
			}

			//check to see if real boundary condition is met (pot at negative elctrode == applied potential)
			electricFieldIntegral = std::accumulate(m_electrostatic[es_electricField].begin(), m_electrostatic[es_electricField].begin() + m_interfacePoint - 5, 0.0) * m_dxf;
			electricFieldIntegral += std::accumulate(m_electrostatic[es_electricField].begin() + m_interfacePoint - 5, m_electrostatic[es_electricField].begin() + m_interfacePoint + m_Isize, 0.0) * m_dxs1;
			electricFieldIntegral += std::accumulate(m_electrostatic[es_electricField].begin() + m_interfacePoint + m_Isize, m_electrostatic[es_electricField].end(), 0.0) * m_dxs2;
			differenceEF = (m_appliedBias - electricFieldIntegral) / m_appliedBias;
			//and apply the updated boundary condition for next round
			electricFieldStart = electricFieldStart * (differenceEF + 1);
		}

		for (array_type::size_type i{ 0 }; i < m_interfacePoint - 4; ++i) // update potential profile using the aquired electrid field
		{
			m_electrostatic[es_potential][i + 1] = m_electrostatic[es_potential][i] - m_electrostatic[es_electricField][i] * m_dxf;
		}


		for (array_type::size_type i{ m_interfacePoint - 4 }; i < (m_interfacePoint + m_Isize + 1); ++i) // update potential profile using the aquired electrid field
		{
			m_electrostatic[es_potential][i + 1] = m_electrostatic[es_potential][i] - m_electrostatic[es_electricField][i] * m_dxs1;
		}

		for (array_type::size_type i{ m_interfacePoint + m_Isize + 1 }; i < (m_size - 1); ++i) // update potential profile using the aquired electrid field
		{
			m_electrostatic[es_potential][i + 1] = m_electrostatic[es_potential][i] - m_electrostatic[es_electricField][i] * m_dxs2;
		}

		//check to see if real boundary condition is met (pot at reference electrode == 0)
		potAtReference = m_electrostatic[es_potential][m_referencePoint];
		//and update the bias to reflect this for the next round
		//divide by relative position of the refernce point (default 1.5) for fastest conversion
		m_appliedBias += potAtReference / m_referencePositionRelative;
	}
}

void DisOxGM::updateRatesGM()
{
	DOS_array::size_type pos{};
	double n{};
	int hwp{ 200 }; //The halfwaypoint of the RRlookup table, were it switches from forward to backward rates
	for (array_type::size_type i{ 1 }; i < m_interfacePoint; ++i)
	{
		n = m_concentrations[carrier_electrons][i];
		if (n > m_RRlookup[0])
		{
			std::cerr << "Electron concentration exceeded maximum! Adjust the reaction rate table to accomodate higher concentrations.";
			std::exit(1);
		}
		pos = static_cast<DOS_array::size_type>(1 + n / m_RRlookup[0] * hwp); //Convert the electron concentration to the correspoding index
		//this is the left position, pos + 1 is the first index where n is higher than the value we are searching for

		m_reactionRates[rate_forward][i] = m_k * (m_RRlookup[pos] + (hwp * n / m_RRlookup[0] - (pos - 1)) * (m_RRlookup[pos + 1] - m_RRlookup[pos]));
		m_reactionRates[rate_backward][i] = m_k* (m_RRlookup[pos + hwp] + (hwp * n / m_RRlookup[0] - (pos - 1)) * (m_RRlookup[pos + hwp + 1] - m_RRlookup[pos + hwp]));
	}
}

void DisOxGM::updateRatesBV()
{
	double kf{ m_k * exp(-(m_electrostatic[es_potential][0] - m_electrostatic[es_potential][1] - m_negativeElectrodeWF + m_E0) / m_energyConvert) };
	double kb{ m_k  *exp((m_electrostatic[es_potential][0] - m_electrostatic[es_potential][1] - m_negativeElectrodeWF + m_E0) / m_energyConvert) };
	m_reactionRates[rate_forward][1] = (kf < 1 ? kf : 1);
	m_reactionRates[rate_backward][1] = (kb < 1 ? kb : 1);
	for (array_type::size_type i{ 2 }; i < m_interfacePoint; ++i)
	{
		kf = m_k * exp(-(m_electrostatic[es_potential][0] - m_electrostatic[es_potential][i] - m_negativeElectrodeWF + m_E0) / m_energyConvert);
		kb = m_k * exp((m_electrostatic[es_potential][0] - m_electrostatic[es_potential][i] - m_negativeElectrodeWF + m_E0) / m_energyConvert);
		m_reactionRates[rate_forward][i] = (kf < 1 ? kf : 1);
		m_reactionRates[rate_backward][i] = (kb < 1 ? kb:1);
	}
}

void DisOxGM::reactGM()
{
	if (++m_stepCounter == 10) //Rates are updated every 10 steps. Can be adjusted here easily to every 5 or every step. 
	{
		updateRatesGM();
		m_stepCounter = 0;
	}


	for (array_type::size_type i{ 1 }; i < m_interfacePoint + 1; ++i)
	{
		m_contaminantCurrents[switch_X][i] = m_contaminants[contaminant_X][i] * m_reactionRates[contaminant_X][i]
			- m_contaminants[contaminant_Xmin][i] * m_reactionRates[contaminant_Xmin][i]; //swtich_x is the amount that goes from X to Xmin
	}

	for (array_type::size_type i{ 1 }; i < m_interfacePoint + 1; ++i)
	{
		m_concentrations[carrier_electrons][i] -= m_contaminantCurrents[switch_X][i];
	}

	for (array_type::size_type i{ 1 }; i < m_interfacePoint + 1; ++i)
	{
		m_contaminants[contaminant_X][i] -= m_contaminantCurrents[switch_X][i];
	}

	for (array_type::size_type i{ 1 }; i < m_interfacePoint + 1; ++i)
	{
		m_contaminants[contaminant_Xmin][i] += m_contaminantCurrents[switch_X][i];
	}
}

void DisOxGM::reactBV()
{
	if (++m_stepCounter == 10) //Rates are updated every 10 steps. Can be adjusted here easily to every 5 or every step. 
	{
		updateRatesBV();
		m_stepCounter = 0;
	}


	for (array_type::size_type i{ 1 }; i < m_interfacePoint-5; ++i)
	{
		m_contaminantCurrents[switch_X][i] = m_contaminants[contaminant_X][i] * m_concentrations[carrier_electrons][i]*1e-26 * m_reactionRates[contaminant_X][i]
			- m_contaminants[contaminant_Xmin][i] * m_reactionRates[contaminant_Xmin][i]; //swtich_x is the amount that goes from X to Xmin
		if (m_contaminantCurrents[switch_X][i] > m_concentrations[carrier_electrons][i]) //Needed only if a first order reaction model is used to prevent negative electron concentration
			m_contaminantCurrents[switch_X][i] = m_concentrations[carrier_electrons][i];
	}

	for (array_type::size_type i{ 1 }; i < m_interfacePoint-5; ++i)
	{
		m_concentrations[carrier_electrons][i] -= m_contaminantCurrents[switch_X][i];
	}

	for (array_type::size_type i{ 1 }; i < m_interfacePoint-5; ++i)
	{
		m_contaminants[contaminant_X][i] -= m_contaminantCurrents[switch_X][i];
	}

	for (array_type::size_type i{ 1 }; i < m_interfacePoint-5; ++i)
	{
		m_contaminants[contaminant_Xmin][i] += m_contaminantCurrents[switch_X][i];
	}
}

void DisOxGM::calculateCurrents()
{
	reactGM(); //reaction takes precedence over currents (arbitrary but much easier this way)
	//electrons first, only up to the interface, using the current function to fill the current array

	for (array_type::size_type i{ 1 }; i < m_interfacePoint - 3; ++i)
	{
		m_currents[carrier_electrons][i] = negativeCurrente(m_concentrations[carrier_electrons][i - 1], m_concentrations[carrier_electrons][i],
			m_currentConstantElectrons, m_electrostatic[es_electricField][i - 1], m_energyConvertxf);
	}

	for (array_type::size_type i{ m_interfacePoint - 3 }; i < m_interfacePoint + 1; ++i)
	{
		m_currents[carrier_electrons][i] = negativeCurrente(m_concentrations[carrier_electrons][i - 1], m_concentrations[carrier_electrons][i],
			m_currentConstantElectrons2, m_electrostatic[es_electricField][i - 1], m_energyConvertxs1);
	}

	//cations next, needs two different loops for the film and solution sections. Ions cannot enter the electrodes, so we dont need to consider the first and last cell
	for (array_type::size_type i{ 2 }; i < m_interfacePoint - 3; ++i)
	{
		m_currents[carrier_cations][i] = positiveCurrent(m_concentrations[carrier_cations][i - 1], m_concentrations[carrier_cations][i],
			m_currentConstantCationsFilm, m_electrostatic[es_electricField][i - 1], m_energyConvertxf);
	}

	for (array_type::size_type i{ m_interfacePoint - 3 }; i < m_interfacePoint + 1; ++i)
	{
		m_currents[carrier_cations][i] = positiveCurrent(m_concentrations[carrier_cations][i - 1], m_concentrations[carrier_cations][i],
			m_currentConstantCationsFilm2, m_electrostatic[es_electricField][i - 1], m_energyConvertxs1);
	}

	for (array_type::size_type i{ m_interfacePoint + 2 }; i < (m_interfacePoint + m_Isize + 1); ++i)
	{
		m_currents[carrier_cations][i] = positiveCurrent(m_concentrations[carrier_cations][i - 1], m_concentrations[carrier_cations][i],
			m_currentConstantCationsFilm2, m_electrostatic[es_electricField][i - 1], m_energyConvertxs1);
	}

	for (array_type::size_type i{ m_interfacePoint + m_Isize + 1 }; i < (m_size - 1); ++i)
	{
		m_currents[carrier_cations][i] = positiveCurrent(m_concentrations[carrier_cations][i - 1], m_concentrations[carrier_cations][i],
			m_currentConstantCationsSolution2, m_electrostatic[es_electricField][i - 1], m_energyConvertxs2);
	}

	//the anions which are very similar to cations
	for (array_type::size_type i{ 2 }; i < m_interfacePoint - 3; ++i)
	{
		m_currents[carrier_anions][i] = negativeCurrent(m_concentrations[carrier_anions][i - 1], m_concentrations[carrier_anions][i],
			m_currentConstantAnionsFilm, m_electrostatic[es_electricField][i - 1], m_energyConvertxf);
	}

	for (array_type::size_type i{ m_interfacePoint - 3 }; i < m_interfacePoint + 1; ++i)
	{

		m_currents[carrier_anions][i] = negativeCurrent(m_concentrations[carrier_anions][i - 1], m_concentrations[carrier_anions][i],
			m_currentConstantAnionsFilm2, m_electrostatic[es_electricField][i - 1], m_energyConvertxs1);
	}

	for (array_type::size_type i{ m_interfacePoint + 2 }; i < (m_interfacePoint + m_Isize + 1); ++i)
	{
		m_currents[carrier_anions][i] = negativeCurrent(m_concentrations[carrier_anions][i - 1], m_concentrations[carrier_anions][i],
			m_currentConstantAnionsFilm2, m_electrostatic[es_electricField][i - 1], m_energyConvertxs1);
	}

	for (array_type::size_type i{ m_interfacePoint + m_Isize + 1 }; i < (m_size - 1); ++i)
	{
		m_currents[carrier_anions][i] = negativeCurrent(m_concentrations[carrier_anions][i - 1], m_concentrations[carrier_anions][i],
			m_currentConstantAnionsSolution2, m_electrostatic[es_electricField][i - 1], m_energyConvertxs2);
	}

	//

	//finally, we need to update the interface currents which we consider the current OVER the interface to be "slow"
	//this means limited by the lower concentration in the film, so we use the solution concentration*QDfillfactor to reflect this
	m_currents[carrier_cations][m_interfacePoint + 1] = positiveCurrent(m_concentrations[carrier_cations][m_interfacePoint], m_concentrations[carrier_cations][m_interfacePoint + 1] * m_QDFillFactor,
		m_currentConstantCationsFilm2, m_electrostatic[es_electricField][m_interfacePoint], m_energyConvertxs1);
	m_currents[carrier_anions][m_interfacePoint + 1] = negativeCurrent(m_concentrations[carrier_anions][m_interfacePoint], m_concentrations[carrier_anions][m_interfacePoint + 1] * m_QDFillFactor,
		m_currentConstantAnionsFilm2, m_electrostatic[es_electricField][m_interfacePoint], m_energyConvertxs1);


	//Added these for X and Xmin
	//Xmin
	
	for (array_type::size_type i{ 2 }; i < m_interfacePoint - 3; ++i)
	{
		m_contaminantCurrents[contaminant_Xmin][i] = negativeCurrent(m_contaminants[contaminant_Xmin][i - 1], m_contaminants[contaminant_Xmin][i],
			m_currentConstantXFilm, m_electrostatic[es_electricField][i - 1], m_energyConvertxf);
	}

	for (array_type::size_type i{ m_interfacePoint - 3 }; i < (m_interfacePoint + m_Isize + 1); ++i)
	{
		m_contaminantCurrents[contaminant_Xmin][i] = negativeCurrent(m_contaminants[contaminant_Xmin][i - 1], m_contaminants[contaminant_Xmin][i],
			m_currentConstantXSolution1, m_electrostatic[es_electricField][i - 1], m_energyConvertxs1);
	}

	for (array_type::size_type i{ m_interfacePoint + m_Isize + 1 }; i < (m_size - 1); ++i)
	{
		m_contaminantCurrents[contaminant_Xmin][i] = negativeCurrent(m_contaminants[contaminant_Xmin][i - 1], m_contaminants[contaminant_Xmin][i],
			m_currentConstantXSolution2, m_electrostatic[es_electricField][i - 1], m_energyConvertxs2);
	}
	//X
	for (array_type::size_type i{ 2 }; i < m_interfacePoint - 3; ++i)
	{
		m_contaminantCurrents[contaminant_X][i] = neutralCurrent(m_contaminants[contaminant_X][i - 1], m_contaminants[contaminant_X][i],
			m_currentConstantXFilm, m_energyConvertxf);
	}

	for (array_type::size_type i{ m_interfacePoint - 3 }; i < (m_interfacePoint + m_Isize + 1); ++i)
	{
		m_contaminantCurrents[contaminant_X][i] = neutralCurrent(m_contaminants[contaminant_X][i - 1], m_contaminants[contaminant_X][i],
			m_currentConstantXSolution1, m_energyConvertxs1);
	}


	for (array_type::size_type i{ m_interfacePoint + m_Isize + 1 }; i < (m_size - 1); ++i)
	{
		m_contaminantCurrents[contaminant_X][i] = neutralCurrent(m_contaminants[contaminant_X][i - 1], m_contaminants[contaminant_X][i],
			m_currentConstantXSolution2, m_energyConvertxs2);
	}
	
	//update total electrons that have entered since the last getCurrent() call:
	m_currentCumulative -= m_currents[carrier_electrons][1];
}

void DisOxGM::updateConcentrations()
{
	//electrons first, only in the film
	for (array_type::size_type i{ 1 }; i < m_interfacePoint + 1; ++i)
	{
		m_concentrations[carrier_electrons][i] += m_currents[carrier_electrons][i] - m_currents[carrier_electrons][i + 1];
	}
	m_concentrations[carrier_electrons][m_interfacePoint - 4] += m_currents[carrier_electrons][m_interfacePoint - 3] * (1 - m_dxs1 / m_dxf);

	//Dirty averaging trick to prevent fast electrons from ruining the simulation
	double ave{ (m_concentrations[carrier_electrons][m_interfacePoint - 3] + m_concentrations[carrier_electrons][m_interfacePoint - 2] + m_concentrations[carrier_electrons][m_interfacePoint - 1]
						+ m_concentrations[carrier_electrons][m_interfacePoint]) / 4 };
	for (array_type::size_type i{ m_interfacePoint - 3 }; i < m_interfacePoint + 1; ++i)
	{
		m_concentrations[carrier_electrons][i] = ave;
	}
	//cations 
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_cations][i] += m_currents[carrier_cations][i] - m_currents[carrier_cations][i + 1];
	}
	m_concentrations[carrier_cations][m_interfacePoint - 4] += m_currents[carrier_cations][m_interfacePoint - 3] * (1 - m_dxs1 / m_dxf);
	m_concentrations[carrier_cations][m_interfacePoint + m_Isize] -= m_currents[carrier_cations][m_interfacePoint + m_Isize] * (1 - m_dxs1 / m_dxs2);
	m_concentrations[carrier_cations][m_referencePoint] = m_saltConcentration;

	//anions
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_anions][i] += m_currents[carrier_anions][i] - m_currents[carrier_anions][i + 1];
	}
	m_concentrations[carrier_anions][m_interfacePoint - 4] += m_currents[carrier_anions][m_interfacePoint - 3] * (1 - m_dxs1 / m_dxf);
	m_concentrations[carrier_anions][m_interfacePoint + m_Isize] -= m_currents[carrier_anions][m_interfacePoint + m_Isize] * (1 - m_dxs1 / m_dxs2);
	m_concentrations[carrier_anions][m_referencePoint] = m_saltConcentration;

	//X
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_contaminants[contaminant_X][i] += m_contaminantCurrents[contaminant_X][i] - m_contaminantCurrents[contaminant_X][i + 1];
	}
	m_contaminants[contaminant_X][m_interfacePoint - 4] += m_contaminantCurrents[contaminant_X][m_interfacePoint - 3] * (1 - m_dxs1 / m_dxf);
	m_contaminants[contaminant_X][m_interfacePoint + m_Isize] -= m_contaminantCurrents[contaminant_X][m_interfacePoint + m_Isize] * (1 - m_dxs1 / m_dxs2);


	//Xmin
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_contaminants[contaminant_Xmin][i] += m_contaminantCurrents[contaminant_Xmin][i] - m_contaminantCurrents[contaminant_Xmin][i + 1];
	}
	m_contaminants[contaminant_Xmin][m_interfacePoint - 4] += m_contaminantCurrents[contaminant_Xmin][m_interfacePoint - 3] * (1 - m_dxs1 / m_dxf);
	m_contaminants[contaminant_Xmin][m_interfacePoint + m_Isize] -= m_contaminantCurrents[contaminant_Xmin][m_interfacePoint + m_Isize] * (1 - m_dxs1 / m_dxs2);

	//X/Xmin equilibriates with the solution at the reference point. Thus an infinite supply of X is available.
	m_contaminants[contaminant_X][m_referencePoint - 1] = m_contaminantConcentration;//m_contaminants[contaminant_Xmin][m_referencePoint - 1];
	m_contaminants[contaminant_Xmin][m_referencePoint - 1] = 0;
}

void DisOxGM::midSave(std::ofstream& midf)
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

	for (int i{ 0 }; i < 2; ++i)
	{
		for (array_type::size_type j{ 0 }; j < m_size; ++j)
		{
			midf << m_reactionRates[i][j] << '\n';
		}
	}

}
