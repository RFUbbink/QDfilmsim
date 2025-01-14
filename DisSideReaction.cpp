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
#include "DisSideReaction.h"

//This is the same class as Discontinuous, but it has additional functionality for calculating the side reactions with molecule 'X', also called contaminant here
//You can set it to use either Butler-Volmer (BV) or Gerischer-Marcus (GM) kinetics, by adjusting the config file. 
//For BV kinetics, you only have to set a k0 (in the config file) and the code will handle the rest.
//For GM kinetics, you do need to provide the ReactionRates file (which you can calculate for certain input parameters using the accompanying python file. 

using array_type = std::vector<std::vector<double>>;

DisSideReaction::DisSideReaction(settings_array& settings, DOS_array& RR)
	:DisCell{ settings },
	m_kineticsType{ (settings[s_kineticsType] < 0.5)? 0 : 1},
	m_contaminantConcentration{ settings[s_contaminantConcentration] },
	m_currentConstantXFilm1{ settings[s_contaminantMobilityFilm] / 1 * settings[s_dt] / m_dxf },
	m_currentConstantXFilm2{ settings[s_contaminantMobilityFilm] / 1 * settings[s_dt] / m_dxs1 },
	m_currentConstantXSolution1{ settings[s_contaminantMobilityFilm] / 1 * settings[s_dt] / m_dxs1 },
	m_currentConstantXSolution2{ settings[s_contaminantMobilitySolution] * settings[s_dt] / m_dxs2 },
	m_k{ m_kineticsType? settings[s_dt]: settings[s_BV_k] * settings[s_dt] },
	m_kIrreversible{(settings[s_kIrriversible] * settings[s_dt] > 1) ? 1: settings[s_kIrriversible] * settings[s_dt] }, // Maxed at 1 to avoid negative concentration of Xmin
	m_E0{ settings[s_BV_E0] },
	m_RRlookup{ RR },
	m_stepCounter{ 0 }

{
	m_contaminants = array_type(3, std::vector<double>(m_size));		//create the oxidadant/reductant concentration array
	m_contaminantCurrents = array_type(3, std::vector<double>(m_size));	//create the reation current array
	m_reactionRates = array_type(3, std::vector<double>(m_size));		//create the array that holds the reaction rates
	m_contaminantsIrreversible = array_type(3, std::vector<double>(m_size));		//create the oxidadant/reductant concentration array
	m_contaminantCurrentsIrreversible = array_type(3, std::vector<double>(m_size));	//create the reation current array
	if (m_kineticsType)
	{
		std::cout << "Currently applying Gerisher-Marcus kinetics\n";
		reactPtr = &DisSideReaction::reactGM;
	}
	else
	{
		std::cout << "Currently applying Butler-Volmer kinetics (or trying anyway)\n";
		reactPtr = &DisSideReaction::reactBV;
	}
	std::cout << "Irriversible kinetics on\n";
}

inline double DisSideReaction::negativeCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//eCon (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (-concentrationLeft * electricField + eCon * (concentrationLeft - concentrationRight)) * curCon;
}

inline double DisSideReaction::negativeCurrente(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon)
{
	//returns the amount of electrons moving to the right cell per m3 per timestep, based on drift/diffusion
	//Has additional functionality that accounts for the electrons occupying higher levels in the DOS of the material, see same function in Discontinuous for explanation. 
	//eCon (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (-concentrationLeft * (electricField - m_electronEnergyFactor * (concentrationLeft - concentrationRight) / pow((concentrationLeft > 0) ? concentrationLeft : 1, 0.3333333333)) + eCon * (concentrationLeft - concentrationRight)) * curCon;
}

inline double DisSideReaction::positiveCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon)
{
	//returns the amount of positively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//eCon (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (concentrationRight * electricField + eCon * (concentrationLeft - concentrationRight)) * curCon;
}

inline double DisSideReaction::neutralCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double eCon)
{
	//returns the amount of neutral particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//eCon (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (eCon * (concentrationLeft - concentrationRight)) * curCon;
}

void DisSideReaction::initializeConcentrations()
{
	//fills the concentration array with cations and anions and contaminant
	for (array_type::size_type i{ 1 }; i < m_interfacePoint + 1; ++i)
	{
		m_concentrations[carrier_cations][i] = m_saltConcentration * m_QDFillFactor;
		m_concentrations[carrier_anions][i] = m_saltConcentration * m_QDFillFactor;
		m_contaminants[contaminant_X][i] = m_contaminantConcentration;
		m_contaminants[contaminant_Xmin][i] = 0; //Setting this to 1 is enough to introduce instability in the case of significant back reaction at V = 0V.
	}

	for (array_type::size_type i{ m_interfacePoint + 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_cations][i] = m_saltConcentration;
		m_concentrations[carrier_anions][i] = m_saltConcentration;
		m_contaminants[contaminant_X][i] = m_contaminantConcentration;
		m_contaminants[contaminant_Xmin][i] = 0;
	}

	m_concentrations[carrier_electrons][0] = 1;
	m_concentrations[carrier_electrons][1] = 1;
}

void DisSideReaction::calculatePotentialProfile()
{
	//This function takes in the reference to ePotentialArray, where the electrostatic potential profile and its 2 space derivates reside
	//It updates the array directly based on the input and thus does not need to return any values.
	m_newAppliedBias = (2 * m_appliedBias - m_oldAppliedBias); // We now predict what the coming applied bias will be based on linear extrapolation from the last 2. 
	//This better estimate allows us to get one-cycle conversion on the applied bias every single time. 
	
	//calculate the second space derivative of the electrostatic potential for every cell based on the Poisson equation
	for (array_type::size_type i{ 1 }; i < m_interfacePoint + 1; ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantFilm * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i]
			- m_concentrations[carrier_anions][i] - m_contaminants[contaminant_Xmin][i] - m_contaminantsIrreversible[contaminant_Xmin][i]);
	}

	for (array_type::size_type i{ m_interfacePoint + 1 }; i < (m_size - 1); ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantSolution * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i]
			- m_concentrations[carrier_anions][i] - m_contaminants[contaminant_Xmin][i] - m_contaminantsIrreversible[contaminant_Xmin][i]);
	}

	double potAtReference{ 1.0 };
	while (std::abs(potAtReference) > 0.0001)
	{
		//electricFieldStart is the initial guess for a boundary condition because the real boundary conditions are impossible to apply directly (assumes that the potential drops linearly over the cell at the start)
		//This boundary condition will then be updated based on the potential profile that is calculated
		//This process is repeated untal a potential profile is generated that satisfies the real boundary conditions (pot at negative electrode == applied potential && pot at reference electrode == 0)
		double electricFieldStart{ m_newAppliedBias/ m_thickness };
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
			differenceEF = (m_newAppliedBias - electricFieldIntegral) / m_newAppliedBias;
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
		m_newAppliedBias += potAtReference / m_referencePositionRelative;
	}
	m_oldAppliedBias = m_appliedBias;
	m_appliedBias = m_newAppliedBias;
}

void DisSideReaction::inspectPotentialODE(std::ofstream& inspectionFile)
{
	//std::cout << m_dxf << '\n' << m_dxs1 << '\n' << m_dxs2 << '\n';
	//This function takes in the reference to ePotentialArray, where the electrostatic potential profile and its 2 space derivates reside
	//It updates the array directly based on the input and thus does not need to return any values. vBias is also updated live to speed up convergence in the next steps

	inspectionFile << "-----Inspection of potential ODE started-----\n";
	inspectionFile << "Current voltage at WE: " << m_electrostatic[es_potential][0] << '\n';
	inspectionFile << "Current applied bias: " << m_appliedBias << '\n';
	unsigned int totalOuterCycles{ 0 };
	unsigned int totalInnerCycles{ 0 };

	m_newAppliedBias = (2 * m_appliedBias - m_oldAppliedBias);
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
	while (std::abs(potAtReference) > 0.0001)
	{
		inspectionFile << "--Started outer cycle--\n";
		totalOuterCycles++;
		//electricFieldStart is the initial guess for a boundary condition because the real boundary conditions are impossible to apply directly (assumes that the potential drops linearly over the cell at the start)
		//This boundary condition will then be updated based on the potential profile that is calculated
		//This process is repeated untal a potential profile is generated that satisfies the real boundary conditions (pot at negative electrode == applied potential && pot at reference electrode == 0)
		double electricFieldStart{ m_newAppliedBias / m_thickness };
		double differenceEF{ 1.0 };
		double electricFieldIntegral{};
		while (std::abs(differenceEF) > 0.00001) //typically takes 2 rounds, very good!
		{
			inspectionFile << "--Started inner cycle--\n";
			totalInnerCycles++;
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
			differenceEF = (m_newAppliedBias - electricFieldIntegral) / m_newAppliedBias;
			inspectionFile << "Difference between applied bias and integral Efield: " << differenceEF << '\n';
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
		inspectionFile << "Potential at reference electrode: " << potAtReference << '\n';
		//and update the bias to reflect this for the next round
		//divide by relative position of the refernce point (default 1.5) for fastest conversion
		m_newAppliedBias += potAtReference / m_referencePositionRelative;
	
	}
	m_oldAppliedBias = m_appliedBias;
	m_appliedBias = m_newAppliedBias;
	inspectionFile << "Outer cycles: " << totalOuterCycles << '\n';
	inspectionFile << "Inner cycles: " << totalInnerCycles << '\n';
	inspectionFile << "-----End of this inspection-----\n\n\n";
}

void DisSideReaction::updateRatesGM()
{
	//Calculates the reaction rates of electrons with molecule X according to Gerischer theory. 
	//The reactionRates lookup table (m_RRlookup) contains the precalculated rates for a given set of input parameters (Calculated in python) for a given concentration of electrons. 
	//Linear interpolation bewteen 2 values in the table to get the rate for a given electron concentration.
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
		m_reactionRates[rate_backward][i] = m_k * (m_RRlookup[pos + hwp] + (hwp * n / m_RRlookup[0] - (pos - 1)) * (m_RRlookup[pos + hwp + 1] - m_RRlookup[pos + hwp]));
	}
}

void DisSideReaction::updateRatesBV()
{
	//Calculates the reaction rates of electrons and molecule X using Butler-Volmer theory
	double kf{ m_k * exp(-(m_electrostatic[es_potential][0] - m_electrostatic[es_potential][1] - m_negativeElectrodeWF + m_E0) / m_energyConvert) };
	double kb{ m_k * exp((m_electrostatic[es_potential][0] - m_electrostatic[es_potential][1] - m_negativeElectrodeWF + m_E0) / m_energyConvert) };
	m_reactionRates[rate_forward][1] = (kf < 1 ? kf : 1);
	m_reactionRates[rate_backward][1] = (kb < 1 ? kb : 1);
	for (array_type::size_type i{ 2 }; i < m_interfacePoint; ++i)
	{
		kf = m_k * exp(-(m_electrostatic[es_potential][0] - m_electrostatic[es_potential][i] - m_negativeElectrodeWF + m_E0) / m_energyConvert);
		kb = m_k * exp((m_electrostatic[es_potential][0] - m_electrostatic[es_potential][i] - m_negativeElectrodeWF + m_E0) / m_energyConvert);
		m_reactionRates[rate_forward][i] = (kf < 1 ? kf : 1);
		m_reactionRates[rate_backward][i] = (kb < 1 ? kb : 1);
	}
}

void DisSideReaction::reactGM()
{
	if (++m_stepCounter == 10) //Rates are updated every 10 steps. Can be adjusted here easily to every 5 or every step. 
	{
		updateRatesGM();
		m_stepCounter = 0;
	}

	//Does the calculation for how much electrons acutally react with X to form X- (or the other way around), taking into account here the concentrations of X and X-
	//Electron concentration was taken into account in the calculation of the rates already
	for (array_type::size_type i{ 1 }; i < m_interfacePoint + 1; ++i)
	{
		m_contaminantCurrents[switch_X][i] = m_contaminants[contaminant_X][i] * m_reactionRates[contaminant_X][i]
			- m_contaminants[contaminant_Xmin][i] * m_reactionRates[contaminant_Xmin][i]; //switch_x is the amount that goes from X to Xmin
		//if (m_contaminantCurrents[switch_X][i] < 0) //Enable these lines to disable the back reaction (Go to ideal irriversible GM behaviour)
		//	m_contaminantCurrents[switch_X][i] = 0;
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

void DisSideReaction::reactBV()
{
	if (++m_stepCounter == 10) //Rates are updated every 10 steps. Can be adjusted here easily to every 5 or every step. 
	{
		updateRatesBV();
		m_stepCounter = 0;
	}

	//Does the calculation for how much electrons acutally react with X to form X- (or the other way around), taking into account here the concentrations of X and X- and also electrons
	for (array_type::size_type i{ 1 }; i < m_interfacePoint - 5; ++i)
	{
		m_contaminantCurrents[switch_X][i] = m_contaminants[contaminant_X][i] * m_concentrations[carrier_electrons][i] * 1e-26 * m_reactionRates[contaminant_X][i]
			- m_contaminants[contaminant_Xmin][i] * m_reactionRates[contaminant_Xmin][i]; //swtich_x is the amount that goes from X to Xmin
		//if (m_contaminantCurrents[switch_X][i] > m_concentrations[carrier_electrons][i]) //Needed only if a first order reaction model is used to prevent negative electron concentration
		//	m_contaminantCurrents[switch_X][i] = m_concentrations[carrier_electrons][i];
	}

	for (array_type::size_type i{ 1 }; i < m_interfacePoint - 5; ++i)
	{
		m_concentrations[carrier_electrons][i] -= m_contaminantCurrents[switch_X][i];
	}

	for (array_type::size_type i{ 1 }; i < m_interfacePoint - 5; ++i)
	{
		m_contaminants[contaminant_X][i] -= m_contaminantCurrents[switch_X][i];
	}

	for (array_type::size_type i{ 1 }; i < m_interfacePoint - 5; ++i)
	{
		m_contaminants[contaminant_Xmin][i] += m_contaminantCurrents[switch_X][i];
	}
}

void DisSideReaction::reactIrreversible()//In case you want irreverisble second step of the X-.
{
	for (array_type::size_type i{ 1 }; i < m_size + 1; ++i)
	{
		m_contaminantCurrentsIrreversible[switch_X][i] = m_contaminants[contaminant_Xmin][i] * m_kIrreversible;
	}

	for (array_type::size_type i{ 1 }; i < m_size + 1; ++i)
	{
		m_contaminants[contaminant_Xmin][i] -= m_contaminantCurrentsIrreversible[switch_X][i];
	}

	for (array_type::size_type i{ 1 }; i < m_size + 1; ++i)
	{
		m_contaminantsIrreversible[contaminant_Xmin][i] += m_contaminantCurrentsIrreversible[switch_X][i];
	}

	//Im also directly handling the movement of the irreverible species here, since it can then be easily turned off by not calling this function to turn off the irriversible behaviour
	//making this species immobile is a bad idea as a LOT of negative charge will build up at the film/electrolyte interface.
	//Properties of this species are copied from regular Xmin: maximum stability of the simulation hopefully and it is not really that important. 

	for (array_type::size_type i{ 2 }; i < m_interfacePoint - 3; ++i)
	{
		m_contaminantCurrentsIrreversible[contaminant_Xmin][i] = negativeCurrent(m_contaminantsIrreversible[contaminant_Xmin][i - 1], m_contaminantsIrreversible[contaminant_Xmin][i],
			m_currentConstantXFilm1, m_electrostatic[es_electricField][i - 1], m_energyConvertxf);
	}

	for (array_type::size_type i{ m_interfacePoint - 3 }; i < m_interfacePoint + 1; ++i)
	{
		m_contaminantCurrentsIrreversible[contaminant_Xmin][i] = negativeCurrent(m_contaminantsIrreversible[contaminant_Xmin][i - 1], m_contaminantsIrreversible[contaminant_Xmin][i],
			m_currentConstantXFilm2, m_electrostatic[es_electricField][i - 1], m_energyConvertxf);
	}

	for (array_type::size_type i{ m_interfacePoint + 1 }; i < (m_interfacePoint + m_Isize + 1); ++i)
	{
		m_contaminantCurrentsIrreversible[contaminant_Xmin][i] = negativeCurrent(m_contaminantsIrreversible[contaminant_Xmin][i - 1], m_contaminantsIrreversible[contaminant_Xmin][i],
			m_currentConstantXSolution1, m_electrostatic[es_electricField][i - 1], m_energyConvertxs1);
	}

	for (array_type::size_type i{ m_interfacePoint + m_Isize + 1 }; i < (m_size - 1); ++i)
	{
		m_contaminantCurrentsIrreversible[contaminant_Xmin][i] = negativeCurrent(m_contaminantsIrreversible[contaminant_Xmin][i - 1], m_contaminantsIrreversible[contaminant_Xmin][i],
			m_currentConstantXSolution2, m_electrostatic[es_electricField][i - 1], m_energyConvertxs2);
	}

	//Xmin
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_contaminantsIrreversible[contaminant_Xmin][i] += m_contaminantCurrentsIrreversible[contaminant_Xmin][i] - m_contaminantCurrentsIrreversible[contaminant_Xmin][i + 1];
	}
	m_contaminantsIrreversible[contaminant_Xmin][m_interfacePoint - 4] += m_contaminantCurrentsIrreversible[contaminant_Xmin][m_interfacePoint - 3] * (1 - m_dxs1 / m_dxf);
	m_contaminantsIrreversible[contaminant_Xmin][m_interfacePoint + m_Isize] -= m_contaminantCurrentsIrreversible[contaminant_Xmin][m_interfacePoint + m_Isize] * (1 - m_dxs1 / m_dxs2);

}

void DisSideReaction::calculateCurrents()
{
	(this->*reactPtr)(); //reaction takes precedence over currents (arbitrary but much easier this way)
	//Also this ugly as sin syntax is the most efficient way to call the proper function. I think. Although if statement may be better?
	//If it was not clear, what happens here is that reactBV or reactGM is called depending on what the user selected in the Config file 
	reactIrreversible();//Then the irreversible reaction happens. Enable or disable based on what you need. 
	
						
	//Then the current densities are calculated below for every carrier, including X and X-
	//Electrons first, up to the film/solution interface 
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


	//Added these for X and X-
	//X-

	for (array_type::size_type i{ 2 }; i < m_interfacePoint - 3; ++i)
	{
		m_contaminantCurrents[contaminant_Xmin][i] = negativeCurrent(m_contaminants[contaminant_Xmin][i - 1], m_contaminants[contaminant_Xmin][i],
			m_currentConstantXFilm1, m_electrostatic[es_electricField][i - 1], m_energyConvertxf);
	}

	for (array_type::size_type i{ m_interfacePoint - 3 }; i < m_interfacePoint + 1; ++i)
	{
		m_contaminantCurrents[contaminant_Xmin][i] = negativeCurrent(m_contaminants[contaminant_Xmin][i - 1], m_contaminants[contaminant_Xmin][i],
			m_currentConstantXFilm2, m_electrostatic[es_electricField][i - 1], m_energyConvertxf);
	}

	for (array_type::size_type i{ m_interfacePoint + 1 }; i < (m_interfacePoint + m_Isize + 1); ++i)
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
			m_currentConstantXFilm1, m_energyConvertxf);
	}

	for (array_type::size_type i{ m_interfacePoint - 3 }; i < m_interfacePoint + 1; ++i)
	{
		m_contaminantCurrents[contaminant_X][i] = neutralCurrent(m_contaminants[contaminant_X][i - 1], m_contaminants[contaminant_X][i],
			m_currentConstantXFilm2, m_energyConvertxf);
	}

	for (array_type::size_type i{ m_interfacePoint + 1 }; i < (m_interfacePoint + m_Isize + 1); ++i)
	{
		m_contaminantCurrents[contaminant_X][i] = neutralCurrent(m_contaminants[contaminant_X][i - 1], m_contaminants[contaminant_X][i],
			m_currentConstantXSolution1, m_energyConvertxs1);
	}


	for (array_type::size_type i{ m_interfacePoint + m_Isize + 1 }; i < (m_size - 1); ++i)
	{
		m_contaminantCurrents[contaminant_X][i] = neutralCurrent(m_contaminants[contaminant_X][i - 1], m_contaminants[contaminant_X][i],
			m_currentConstantXSolution2, m_energyConvertxs2);
	}

	m_currents[carrier_cations][m_interfacePoint + 1] = positiveCurrent(m_concentrations[carrier_cations][m_interfacePoint], m_concentrations[carrier_cations][m_interfacePoint + 1] * m_QDFillFactor,
		m_currentConstantCationsFilm2, m_electrostatic[es_electricField][m_interfacePoint], m_energyConvertxs1);
	m_currents[carrier_anions][m_interfacePoint + 1] = negativeCurrent(m_concentrations[carrier_anions][m_interfacePoint], m_concentrations[carrier_anions][m_interfacePoint + 1] * m_QDFillFactor,
		m_currentConstantAnionsFilm2, m_electrostatic[es_electricField][m_interfacePoint], m_energyConvertxs1);

	//update total electrons that have entered since the last getCurrent() call:
	m_currentCumulative -= m_currents[carrier_electrons][1];
}

void DisSideReaction::updateConcentrations()
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

void DisSideReaction::midSave(std::ofstream& midf)
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
			midf << m_contaminantsIrreversible[i][j] << '\n';
		}
	}

	for (int i{ 0 }; i < 2; ++i)
	{
		for (array_type::size_type j{ 0 }; j < m_size; ++j)
		{
			midf << m_contaminantCurrentsIrreversible[i][j] << '\n';
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
