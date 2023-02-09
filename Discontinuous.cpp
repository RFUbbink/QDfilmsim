#include <cmath>
#include <cstddef>
#include <vector>
#include <array>
#include <numeric>
#include <iostream>
#include <fstream>
#include "Config.h"
#include "Discontinuous.h" 

using array_type = std::vector<std::vector<double>>;

//Member functions of Cell class:

DisCell::DisCell(settings_array& settings) //initialize all the constants necessary for operation, using settings array
//vBias is the TOTAL potential drop over the entire system == the voltage difference between the working and counter electrode
//an estimation of the initial voltage drop over the counterelectrode is made to calculate the initial Vb.
	:m_appliedBias{ settings[s_startVoltage] + settings[s_startVoltage] * settings[s_workingElectrodeArea]
						/ settings[s_counterElectrodeArea] * settings[s_epsilonrFilm] / settings[s_epsilonrSolution] },
	m_voltageIncrement{ settings[s_voltageIncrement] },
	m_saltConcentration{ settings[s_ionConcentration] },
	m_interfacePoint{ static_cast<array_type::size_type>(30) },
	m_referencePoint{ static_cast<array_type::size_type>(m_size-200) },
	m_referencePositionRelative{ (settings[s_cellThickness] + settings[s_filmThickness]) / 2 / settings[s_cellThickness]},
	m_thickness{ settings[s_cellThickness] },
	m_dxf{ (settings[s_filmThickness]-25e-9) / (m_interfacePoint-5) },
	m_dxs1{ 5e-9 }, //adapt the above and below (-etc) values together with this one to change interface resolution
	m_dxs2{ (settings[s_cellThickness] - 300e-9 - settings[s_filmThickness]) / (m_size - m_interfacePoint - m_Isize - 1) },
	m_injectionBarrier{ settings[s_LUMO] - settings[s_negativeElectrodeWF] },
	m_densityOfStates{ settings[s_densityOfStates] },
	m_LUMO{ settings[s_LUMO] },
	m_negativeElectrodeWF{ settings[s_negativeElectrodeWF] },
	m_QDFillFactor{ 1 - settings[s_QDspacefill] },
	m_electronEnergyFactor{ settings[39] / m_dxf },
	m_currentConstantElectrons{ settings[s_electronMobility] * settings[s_dt] / m_dxf },
	m_currentConstantElectrons2{ settings[s_electronMobility] * settings[s_dt] / m_dxs1 },
	m_currentConstantCationsFilm{ settings[s_cationMobilityFilm] * settings[s_dt] / m_dxf },
	m_currentConstantCationsFilm2{ settings[s_cationMobilityFilm] * settings[s_dt] / m_dxs1 },
	//Weird quirk on where the mobility of ions is reduced when they near the interface
	m_currentConstantCationsSolution1{ settings[s_cationMobilityFilm] * settings[s_dt] / m_dxs1 },
	m_currentConstantCationsSolution2{ settings[s_cationMobilitySolution] * settings[s_dt] / m_dxs2 },
	m_currentConstantAnionsFilm{ settings[s_anionMobilityFilm] * settings[s_dt] / m_dxf },
	m_currentConstantAnionsFilm2{ settings[s_anionMobilityFilm] * settings[s_dt] / m_dxs1 },
	m_currentConstantAnionsSolution1{ settings[s_anionMobilityFilm] * settings[s_dt] / m_dxs1 },
	m_currentConstantAnionsSolution2{ settings[s_anionMobilitySolution] * settings[s_dt] / m_dxs2 },
	m_energyConvert{ phys::k * settings[s_temperature] / phys::q },
	m_energyConvertxf{ phys::k * settings[s_temperature] / phys::q / m_dxf },
	m_energyConvertxs1{ phys::k * settings[s_temperature] / phys::q / m_dxs1 },
	m_energyConvertxs2{ phys::k * settings[s_temperature] / phys::q / m_dxs2 },
	m_poissonConstantFilm{ -phys::q / phys::eps0 / settings[s_epsilonrFilm] },
	m_poissonConstantSolution{ -phys::q / phys::eps0 / settings[s_epsilonrSolution] },
	m_currentConvert{ phys::q * m_dxf / settings[s_dt] / 10000 }
{
	m_electrostatic = array_type(3, std::vector<double>(m_size));		//create the potential array
	m_electrostatic[es_potential][0] = settings[s_startVoltage];					//apply the inital bias: updated everytime there is a voltage increment.
																		//This array value should be constant unless there is an increment
	m_concentrations = array_type(3, std::vector<double>(m_size));		//create the concentration array
	m_currents = array_type(3, std::vector<double>(m_size));			//create the current array
	std::cout << "Current interface resolution: " << m_dxs1*1e9 << " nm\n";
}

double DisCell::getCurrent() //returns the electron current at a given time (measured in the first cell next to the negative electrode)
{
	double curcur{ m_currentCumulative * m_currentConvert };
	m_currentCumulative = 0;
	return curcur;
}

double DisCell::getLeakCurrent() //returns the electron current at a given time (measured in the first cell next to the negative electrode)
{
	double curcur{ m_leakCurrentCumulative * m_currentConvert };
	m_leakCurrentCumulative = 0;
	return curcur;
}

double DisCell::getVoltage()
{
	return m_electrostatic[es_potential][0];
}

void DisCell::injectElectrons(const DOS_array& DOS)
{
	//These lines use the real DOS and Boltzmann-distribution to calculate the electron concentration.
	//DOS array also holds the "official" LUMO energy and energy integrating step dE
	double n{};
	double energy{ DOS[0]};
	double dE{ DOS[1] };
	//integrating from Estart to the end of DOS eV using steps of 0.01 eV
	for (DOS_array::size_type i{ 3 }; i < 300; ++i)
	{
		n += DOS[i] / (exp((energy - (m_negativeElectrodeWF - m_electrostatic[es_potential][0] + m_electrostatic[es_potential][1])) / m_energyConvert) + 1) * dE;
		energy += dE;
	}
	//m_concentrations[carrier_electrons][0] = n / 15;
	m_concentrations[carrier_electrons][0] = m_concentrations[carrier_electrons][0] * n / m_concentrations[carrier_electrons][1];
}

void DisCell::calculatePotentialProfile()
{
	//std::cout << m_dxf << '\n' << m_dxs1 << '\n' << m_dxs2 << '\n';
	//This function takes in the reference to ePotentialArray, where the electrostatic potential profile and its 2 space derivates reside
	//It updates the array directly based on the input and thus does not need to return any values. vBias is also updated live to speed up convergence in the next steps

	//calculate the second space derivative of the electrostatic potential for every cell based on the Poisson equation
	for (array_type::size_type i{ 1 }; i < m_interfacePoint + 1; ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantFilm * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i] - m_concentrations[carrier_anions][i]);
	}

	for (array_type::size_type i{ m_interfacePoint + 1}; i < (m_size - 1); ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantSolution * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i] - m_concentrations[carrier_anions][i]);
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
			for (array_type::size_type i{ 0 }; i < m_interfacePoint-5; ++i) //build electric field sequentially using said boundary condition
			{
				m_electrostatic[es_electricField][i + 1] = m_electrostatic[es_electricField][i] - m_electrostatic[es_secondDerivative][i + 1] * m_dxf;
			}


			for (array_type::size_type i{ m_interfacePoint-5}; i < (m_interfacePoint + m_Isize); ++i) //build electric field sequentially using said boundary condition
			{
				m_electrostatic[es_electricField][i + 1] = m_electrostatic[es_electricField][i] - m_electrostatic[es_secondDerivative][i + 1] * m_dxs1;
			}

			for (array_type::size_type i{ m_interfacePoint + m_Isize }; i < (m_size - 2); ++i) //build electric field sequentially using said boundary condition
			{
				m_electrostatic[es_electricField][i + 1] = m_electrostatic[es_electricField][i] - m_electrostatic[es_secondDerivative][i + 1] * m_dxs2;
			}

			//check to see if real boundary condition is met (pot at negative elctrode == applied potential)
			electricFieldIntegral = std::accumulate(m_electrostatic[es_electricField].begin(), m_electrostatic[es_electricField].begin() + m_interfacePoint-5, 0.0) * m_dxf;
			electricFieldIntegral += std::accumulate(m_electrostatic[es_electricField].begin() + m_interfacePoint-5, m_electrostatic[es_electricField].begin() + m_interfacePoint + m_Isize, 0.0) * m_dxs1;
			electricFieldIntegral += std::accumulate(m_electrostatic[es_electricField].begin() + m_interfacePoint + m_Isize, m_electrostatic[es_electricField].end(), 0.0) * m_dxs2;
			differenceEF = (m_appliedBias - electricFieldIntegral) / m_appliedBias;
			//and apply the updated boundary condition for next round
			electricFieldStart = electricFieldStart * (differenceEF + 1);
		}

		for (array_type::size_type i{ 0 }; i < m_interfacePoint-4; ++i) // update potential profile using the aquired electrid field
		{
			m_electrostatic[es_potential][i + 1] = m_electrostatic[es_potential][i] - m_electrostatic[es_electricField][i] * m_dxf;
		}


		for (array_type::size_type i{ m_interfacePoint -4}; i < (m_interfacePoint + m_Isize + 1); ++i) // update potential profile using the aquired electrid field
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

void DisCell::initializeConcentrations()
{
	std::cout << "Electron averaging on!\nSlow mobility near interface on!\n";
	//fills the concentration array with cations and anions
	for (array_type::size_type i{ 1 }; i < m_interfacePoint + 1; ++i)
	{
		m_concentrations[carrier_cations][i] = m_saltConcentration * m_QDFillFactor;
		m_concentrations[carrier_anions][i] = m_saltConcentration* m_QDFillFactor;
	}

	for (array_type::size_type i{ m_interfacePoint + 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_cations][i] = m_saltConcentration;
		m_concentrations[carrier_anions][i] = m_saltConcentration;
	}

	m_concentrations[carrier_electrons][0] = 1;
	m_concentrations[carrier_electrons][1] = 1;
}

inline double DisCell::negativeCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (-concentrationLeft * electricField + eCon * (concentrationLeft - concentrationRight)) * curCon;
}

inline double DisCell::negativeCurrente(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon)
{
	//returns the amount of negatively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	//if (concentrationLeft > 1e25)
	return (-concentrationLeft * (electricField - m_electronEnergyFactor * (concentrationLeft - concentrationRight) / pow((concentrationLeft>0) ? concentrationLeft : 1, 0.3333333333)) + eCon * (concentrationLeft - concentrationRight)) * curCon;
	//else
	//	return (-concentrationLeft * electricField  + m_energyConvertx * (concentrationLeft - concentrationRight)) * curCon;
}

inline double DisCell::positiveCurrent(const double concentrationLeft, const double concentrationRight, const double curCon, const double electricField, const double eCon)
{
	//returns the amount of positively charged particles moving to the right cell per m3 per timestep, based on drift/diffusion
	//energyConvertX (merged for speed) = k*T/q/dx
	//curCon (merged for speed)		= mobility*dt/dx
	return (concentrationRight * electricField + eCon * (concentrationLeft - concentrationRight)) * curCon;
}

void DisCell::calculateCurrents()
{
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
	for (array_type::size_type i{ 1 }; i < m_interfacePoint -3; ++i)
	{
		/*
		if (i == 10)
		{
			m_currents[carrier_electrons][i] = negativeCurrente(m_concentrations[carrier_electrons][i - 1], m_concentrations[carrier_electrons][i],
				m_currentConstantElectrons, m_electrostatic[es_electricField][i - 1]-0.01/m_dxf, m_energyConvertxf);
		}
		else if (i == 15)
		{
			m_currents[carrier_electrons][i] = negativeCurrente(m_concentrations[carrier_electrons][i - 1], m_concentrations[carrier_electrons][i],
				m_currentConstantElectrons, m_electrostatic[es_electricField][i - 1]+0.01/m_dxf, m_energyConvertxf);
		}
		else
		*/
		{
			m_currents[carrier_electrons][i] = negativeCurrente(m_concentrations[carrier_electrons][i - 1], m_concentrations[carrier_electrons][i],
				m_currentConstantElectrons, m_electrostatic[es_electricField][i - 1], m_energyConvertxf);
		}
	}

	//These lines prevent the simulation from crashing at nearly the end, when the extraction takes place. 
	//They prevent the electron concetration from going under zero in near the interface during extraction.
	if (m_currents[carrier_electrons][m_interfacePoint-4] < -m_concentrations[carrier_electrons][m_interfacePoint - 4]) 
		m_currents[carrier_electrons][m_interfacePoint - 4] = -m_concentrations[carrier_electrons][m_interfacePoint - 4];

	for (array_type::size_type i{ m_interfacePoint - 3 }; i < m_interfacePoint +1; ++i)
	{
		m_currents[carrier_electrons][i] = negativeCurrente(m_concentrations[carrier_electrons][i - 1], m_concentrations[carrier_electrons][i],
			m_currentConstantElectrons2, m_electrostatic[es_electricField][i - 1], m_energyConvertxs1);
		//The following if statement makes sure that the concentraton can never be negative by limiting the current density in problematic areas
		//Current is then reduced (both ways, since both can be problematic) This is a bit of a "timeskip" trick, not entirely physical, but neither are negative concentrations
		//at lower dt, this should never have an impact
		if (m_currents[carrier_electrons][i] - m_currents[carrier_electrons][i - 1] > m_concentrations[carrier_electrons][i])
		{
			m_currents[carrier_electrons][i-1] = m_currents[carrier_electrons][i-1] / (m_currents[carrier_electrons][i] - m_currents[carrier_electrons][i - 1]) * m_concentrations[carrier_electrons][i];
			m_currents[carrier_electrons][i] = m_currents[carrier_electrons][i] / (m_currents[carrier_electrons][i] - m_currents[carrier_electrons][i - 1]) * m_concentrations[carrier_electrons][i];
		}
	}

	//double cons{};
	//double mod{};
	//cations next, needs two different loops for the film and solution sections. Ions cannot enter the electrodes, so we dont need to consider the first and last cell
	for (array_type::size_type i{ 2 }; i < m_interfacePoint - 3; ++i)
	{
		//cons = sqrt(m_concentrations[carrier_cations][i]/6e26 + 1e-10);
		//mod = exp(-30 * ((cons) / (1 + cons) - 0.3 * m_concentrations[carrier_cations][i]/6e26));
		m_currents[carrier_cations][i] = positiveCurrent(m_concentrations[carrier_cations][i - 1], m_concentrations[carrier_cations][i],
			m_currentConstantCationsFilm, m_electrostatic[es_electricField][i - 1], m_energyConvertxf);
	}

	for (array_type::size_type i{ m_interfacePoint - 3 }; i < m_interfacePoint + 1; ++i)
	{
		//cons = sqrt(m_concentrations[carrier_cations][i]/6e26 + 1e-10);
		//mod = exp(-30 * ((cons) / (1 + cons) - 0.3 * m_concentrations[carrier_cations][i]/6e26));
		m_currents[carrier_cations][i] = positiveCurrent(m_concentrations[carrier_cations][i - 1], m_concentrations[carrier_cations][i],
			m_currentConstantCationsFilm2, m_electrostatic[es_electricField][i - 1], m_energyConvertxs1);
	}

	for (array_type::size_type i{ m_interfacePoint + 2 }; i < (m_interfacePoint + m_Isize + 1); ++i)
	{
		m_currents[carrier_cations][i] = positiveCurrent(m_concentrations[carrier_cations][i - 1], m_concentrations[carrier_cations][i],
			m_currentConstantCationsSolution1, m_electrostatic[es_electricField][i - 1], m_energyConvertxs1);
	}

	for (array_type::size_type i{ m_interfacePoint + m_Isize + 1 }; i < (m_size - 1); ++i)
	{
		m_currents[carrier_cations][i] = positiveCurrent(m_concentrations[carrier_cations][i - 1], m_concentrations[carrier_cations][i],
			m_currentConstantCationsSolution2, m_electrostatic[es_electricField][i - 1], m_energyConvertxs2);
	}

	//the anions which are very similar to cations
	for (array_type::size_type i{ 2 }; i < m_interfacePoint -3; ++i)
	{
		//cons = sqrt(m_concentrations[carrier_anions][i]/6e26+1e-10);
		//mod = exp(-30 * ((cons) / (1 + cons) - 0.3 * m_concentrations[carrier_anions][i]/6e26));
		m_currents[carrier_anions][i] = negativeCurrent(m_concentrations[carrier_anions][i - 1], m_concentrations[carrier_anions][i],
			m_currentConstantAnionsFilm, m_electrostatic[es_electricField][i - 1], m_energyConvertxf);
	}

	for (array_type::size_type i{ m_interfacePoint-3 }; i < m_interfacePoint + 1; ++i)
	{
		//cons = sqrt(m_concentrations[carrier_anions][i]/6e26+1e-10);
		//mod = exp(-30 * ((cons) / (1 + cons) - 0.3 * m_concentrations[carrier_anions][i]/6e26));
		m_currents[carrier_anions][i] = negativeCurrent(m_concentrations[carrier_anions][i - 1], m_concentrations[carrier_anions][i],
			m_currentConstantAnionsFilm2, m_electrostatic[es_electricField][i - 1], m_energyConvertxs1);
	}

	for (array_type::size_type i{ m_interfacePoint + 2 }; i < (m_interfacePoint + m_Isize + 1); ++i)
	{
		m_currents[carrier_anions][i] = negativeCurrent(m_concentrations[carrier_anions][i - 1], m_concentrations[carrier_anions][i],
			m_currentConstantAnionsSolution1, m_electrostatic[es_electricField][i - 1], m_energyConvertxs1);
	}

	for (array_type::size_type i{ m_interfacePoint + m_Isize + 1 }; i < (m_size - 1); ++i)
	{
		m_currents[carrier_anions][i] = negativeCurrent(m_concentrations[carrier_anions][i - 1], m_concentrations[carrier_anions][i],
			m_currentConstantAnionsSolution2, m_electrostatic[es_electricField][i - 1], m_energyConvertxs2);
	}

	//finally, we need to update the interface currents which we consider the current OVER the interface to be "slow"
	//this means limited by the lower concentration in the film, so we use the solution concentration*QDfillfactor to reflect this

	//SQRT implementation
	//double catcur{ (sqrt(m_concentrations[carrier_cations][m_interfacePoint]* m_concentrations[carrier_cations][m_interfacePoint+1]* m_QDFillFactor) * m_electrostatic[es_electricField][m_interfacePoint] +
	//	m_energyConvertxf * (m_concentrations[carrier_cations][m_interfacePoint] - m_concentrations[carrier_cations][m_interfacePoint+1] * m_QDFillFactor)) * m_currentConstantCationsFilm };
	//m_currents[carrier_cations][m_interfacePoint+1] = (catcur < m_concentrations[carrier_cations][m_interfacePoint+1] ? catcur : m_concentrations[carrier_cations][m_interfacePoint+1]);
	//if (m_currents[carrier_cations][m_interfacePoint] < (m_concentrations[carrier_cations][m_interfacePoint - 1] - 2e26))
	//	m_currents[carrier_cations][m_interfacePoint] = (m_concentrations[carrier_cations][m_interfacePoint - 1] - 2e26);

	//double ancur{ (-sqrt(m_concentrations[carrier_anions][m_interfacePoint] * m_concentrations[carrier_anions][m_interfacePoint + 1] * m_QDFillFactor) * m_electrostatic[es_electricField][m_interfacePoint] +
	//	m_energyConvertxf * (m_concentrations[carrier_anions][m_interfacePoint] - m_concentrations[carrier_anions][m_interfacePoint+1] * m_QDFillFactor)) * m_currentConstantAnionsFilm };
    //m_currents[carrier_anions][m_interfacePoint + 1] = (ancur < m_concentrations[carrier_anions][m_interfacePoint+1] ? ancur : m_concentrations[carrier_anions][m_interfacePoint+1]);

	//BOLTZMANN implementation
	//double catcur{ (m_concentrations[carrier_cations][m_interfacePoint] - m_concentrations[carrier_cations][m_interfacePoint+1] * m_QDFillFactor *
	//												exp((m_electrostatic[es_potential][m_interfacePoint+1] - m_electrostatic[es_potential][m_interfacePoint]) / m_energyConvert))/10000 };
	//m_currents[carrier_cations][m_interfacePoint+1] = (catcur < m_concentrations[carrier_cations][m_interfacePoint+1] ? catcur : m_concentrations[carrier_cations][m_interfacePoint+1]);
	//double ancur{(m_concentrations[carrier_anions][m_interfacePoint] - m_concentrations[carrier_anions][m_interfacePoint+1] * m_QDFillFactor *
	//												exp((m_electrostatic[es_potential][m_interfacePoint] - m_electrostatic[es_potential][m_interfacePoint+1]) / m_energyConvert))/10000 };
	//m_currents[carrier_anions][m_interfacePoint+1] = (ancur < m_concentrations[carrier_anions][m_interfacePoint+1] ? ancur : m_concentrations[carrier_anions][m_interfacePoint+1]);


	//double catcur{ (m_concentrations[carrier_cations][m_interfacePoint] * m_QDFillFactor * 3 * m_electrostatic[es_electricField][m_interfacePoint - 1] +
	//	m_energyConvertx * (m_concentrations[carrier_cations][m_interfacePoint - 1] - m_concentrations[carrier_cations][m_interfacePoint] * m_QDFillFactor)) * m_currentConstantCationsFilm };
	//m_currents[carrier_cations][m_interfacePoint] = (catcur < m_concentrations[carrier_cations][m_interfacePoint] ? catcur : m_concentrations[carrier_cations][m_interfacePoint]);

	m_currents[carrier_cations][m_interfacePoint + 1] = positiveCurrent(m_concentrations[carrier_cations][m_interfacePoint], m_concentrations[carrier_cations][m_interfacePoint + 1] * m_QDFillFactor,
		m_currentConstantCationsFilm2, m_electrostatic[es_electricField][m_interfacePoint], m_energyConvertxs1);
	m_currents[carrier_anions][m_interfacePoint+1] = negativeCurrent(m_concentrations[carrier_anions][m_interfacePoint], m_concentrations[carrier_anions][m_interfacePoint + 1] * m_QDFillFactor,
		m_currentConstantAnionsFilm2, m_electrostatic[es_electricField][m_interfacePoint], m_energyConvertxs1);


	//update total electrons that have entered since the last getCurrent() call:
	m_currentCumulative -= m_currents[carrier_electrons][1];

}

void DisCell::updateConcentrations()
{
	//electrons first, only in the film
	for (array_type::size_type i{ 1 }; i < m_interfacePoint + 1; ++i)
	{
		m_concentrations[carrier_electrons][i] += m_currents[carrier_electrons][i] - m_currents[carrier_electrons][i + 1];
	}
	m_concentrations[carrier_electrons][m_interfacePoint - 4] += m_currents[carrier_electrons][m_interfacePoint - 3] * (1 - m_dxs1 / m_dxf);

	//Dirty averaging trick to prevent fast electrons from ruining the simulation
	double ave{ (m_concentrations[carrier_electrons][m_interfacePoint - 3] + m_concentrations[carrier_electrons][m_interfacePoint - 2] + m_concentrations[carrier_electrons][m_interfacePoint - 1]
						+ m_concentrations[carrier_electrons][m_interfacePoint])/4};
	for (array_type::size_type i{ m_interfacePoint - 3 }; i < m_interfacePoint + 1; ++i)
	{
		m_concentrations[carrier_electrons][i] = ave;
	}
	//cations 
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_cations][i] += m_currents[carrier_cations][i] - m_currents[carrier_cations][i + 1];
	}
	m_concentrations[carrier_cations][m_interfacePoint-4] += m_currents[carrier_cations][m_interfacePoint -3] * (1 - m_dxs1 / m_dxf);
	m_concentrations[carrier_cations][m_interfacePoint + m_Isize] -= m_currents[carrier_cations][m_interfacePoint + m_Isize] * (1 - m_dxs1 / m_dxs2);
	m_concentrations[carrier_cations][m_referencePoint] = m_saltConcentration;

	//anions
	for (array_type::size_type i{ 1 }; i < (m_size - 1); ++i)
	{
		m_concentrations[carrier_anions][i] += m_currents[carrier_anions][i] - m_currents[carrier_anions][i + 1];
	}
	m_concentrations[carrier_anions][m_interfacePoint-4] += m_currents[carrier_anions][m_interfacePoint -3] * (1 - m_dxs1 / m_dxf);
	m_concentrations[carrier_anions][m_interfacePoint + m_Isize] -= m_currents[carrier_anions][m_interfacePoint + m_Isize] * (1 - m_dxs1 / m_dxs2);
	m_concentrations[carrier_anions][m_referencePoint] = m_saltConcentration;

}

void DisCell::changeBias(double vBiasChange)
{
	m_appliedBias += vBiasChange;
	m_electrostatic[es_potential][0] += vBiasChange;
} //increases the applied bias by vbuaschange

void DisCell::resetInjection()
{
	if (m_concentrations[carrier_electrons][0] < 1)
		m_concentrations[carrier_electrons][0] = 1;
	if (m_concentrations[carrier_electrons][1] < 1)
		m_concentrations[carrier_electrons][1] = 1;
}

void DisCell::adjustMobility(settings_array& settings)
{
	if (m_concentrations[carrier_cations][2] > 1e26)
	{
		m_currentConstantCationsFilm = settings[s_cationMobilityFilm] * settings[s_dt] / m_dxf;
		m_currentConstantAnionsFilm = settings[s_anionMobilityFilm] * settings[s_dt] / m_dxf;
	}
	m_currentConstantCationsFilm = (settings[s_cationMobilityFilm] + (1 - m_concentrations[carrier_cations][2]/1e26)*(settings[s_cationMobilitySolution]- settings[s_cationMobilityFilm])) * settings[s_dt] / m_dxf;
	m_currentConstantAnionsFilm = (settings[s_anionMobilityFilm] + (1 - m_concentrations[carrier_cations][2] / 1e26)* (settings[s_anionMobilitySolution] - settings[s_anionMobilityFilm])) * settings[s_dt] / m_dxf;
}

DisCell& DisCell::operator++()
{
	//manage the applied bias by incrementing (potentialArray is the big one, appliedBias is also incrementen just to speed up convergense in calculatePotentialProfile)
	m_appliedBias += m_voltageIncrement;
	m_electrostatic[es_potential][0] += m_voltageIncrement;
	return *this;
}

DisCell& DisCell::operator--()
{
	//manage the applied bias by incrementing (potentialArray is the big one, appliedBias is also incrementen just to speed up convergense in calculatePotentialProfile)
	m_appliedBias -= m_voltageIncrement;
	m_electrostatic[es_potential][0] -= m_voltageIncrement;
	return *this;
}

void DisCell::midSave(std::ofstream& midf)
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
