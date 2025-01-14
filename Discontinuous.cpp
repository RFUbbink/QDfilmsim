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

//This is the same class as Electrochemistry, but uses discontinuous space intervalls. 
//Many functions are the same, so if not clear here there could be comments that clatify them in Electrochemistry.cpp

DisCell::DisCell(settings_array& settings) //initialize all the constants necessary for operation, using settings array
//vBias is the TOTAL potential drop over the entire system == the voltage difference between the working and counter electrode
//an estimation of the initial voltage drop over the counterelectrode is made to calculate the initial Vb.
	:m_appliedBias{ settings[s_startVoltage] + settings[s_startVoltage] * settings[s_epsilonrFilm] / settings[s_epsilonrSolution] },
	m_oldAppliedBias{ settings[s_startVoltage] + settings[s_startVoltage] * settings[s_epsilonrFilm] / settings[s_epsilonrSolution] },
	m_voltageIncrement{ settings[s_voltageIncrement] },
	m_saltConcentration{ settings[s_ionConcentration] },
	m_size{ static_cast<array_type::size_type>(settings[s_amountOfCells]) },
	m_Isize{ static_cast<array_type::size_type>(settings[s_amountofInterfaceCells]) },
	m_interfacePoint{ static_cast<array_type::size_type>(settings[s_amountOfFilmCells]) },
	m_referencePoint{ static_cast<array_type::size_type>(m_size - (m_size-m_interfacePoint-m_Isize)/2) }, //So reference distance in the config does nothing here
	m_thickness{ settings[s_cellThickness] },
	m_dxs1{ settings[s_interfaceResolution] * 1e-9},
	m_dxf{ (settings[s_filmThickness] - 5 * m_dxs1) / (m_interfacePoint-5) }, 
	m_dxs2{ (settings[s_cellThickness] - m_Isize * m_dxs1 - settings[s_filmThickness]) / (m_size - m_interfacePoint - m_Isize - 1) },
	m_referencePositionRelative{ (settings[s_cellThickness] - ((m_size - m_interfacePoint - m_Isize) / 2 - 1) * m_dxs2 ) / settings[s_cellThickness] },
	m_injectionBarrier{ settings[s_LUMO] - settings[s_negativeElectrodeWF] },
	m_LUMO{ settings[s_LUMO] },
	m_negativeElectrodeWF{ settings[s_negativeElectrodeWF] },
	m_QDFillFactor{ 1 - settings[s_QDspacefill] },
	m_electronEnergyFactor{ settings[39] / m_dxf },
	m_currentConstantElectrons{ settings[s_electronMobility] * settings[s_dt] / m_dxf },
	m_currentConstantElectrons2{ settings[s_electronMobility] * settings[s_dt] / m_dxs1 },
	m_currentConstantCationsFilm{ settings[s_cationMobilityFilm] * settings[s_dt] / m_dxf },
	m_currentConstantCationsFilm2{ settings[s_cationMobilityFilm] * settings[s_dt] / m_dxs1 },
	//Weird quirk on where the mobility of ions is reduced when they near the interface (in the interface region)
	m_currentConstantCationsSolution1{ settings[s_cationMobilityFilm] * settings[s_dt] / m_dxs1 },
	m_currentConstantCationsSolution2{ settings[s_cationMobilitySolution] * settings[s_dt] / m_dxs2 },
	m_currentConstantAnionsFilm{ settings[s_anionMobilityFilm] * settings[s_dt] / m_dxf },
	m_currentConstantAnionsFilm2{ settings[s_anionMobilityFilm] * settings[s_dt] / m_dxs1 },
	m_currentConstantAnionsSolution1{ settings[s_anionMobilityFilm] * settings[s_dt] / m_dxs1 },
	m_currentConstantAnionsSolution2{ settings[s_anionMobilitySolution] * settings[s_dt] / m_dxs2 },
	m_energyConvert{ physics::k * settings[s_temperature] / physics::q },
	m_energyConvertxf{ physics::k * settings[s_temperature] / physics::q / m_dxf },
	m_energyConvertxs1{ physics::k * settings[s_temperature] / physics::q / m_dxs1 },
	m_energyConvertxs2{ physics::k * settings[s_temperature] / physics::q / m_dxs2 },
	m_poissonConstantFilm{ -physics::q / physics::eps0 / settings[s_epsilonrFilm] },
	m_poissonConstantSolution{ -physics::q / physics::eps0 / settings[s_epsilonrSolution] },
	m_currentConvert{ physics::q * m_dxf / settings[s_dt] / 10000 }
{
	m_electrostatic = array_type(3, std::vector<double>(m_size));		//create the potential array
	m_electrostatic[es_potential][0] = settings[s_startVoltage];					//apply the inital bias: updated everytime there is a voltage increment.
																		//This array value should be constant unless there is an increment
	m_concentrations = array_type(3, std::vector<double>(m_size));		//create the concentration array
	m_currents = array_type(3, std::vector<double>(m_size));			//create the current array
	std::cout << "Current interface resolution: " << m_dxs1*1e9 << " nm\n";
	std::cout << "Electron averaging on!\nSlow mobility near interface on!\n";
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
	//This function takes in the reference to ePotentialArray, where the electrostatic potential profile and its 2 space derivates reside
	//It updates the array directly based on the input and thus does not need to return any values. m_appliedBias is also updated live to speed up convergence in the next steps
	m_newAppliedBias = (2 * m_appliedBias - m_oldAppliedBias); //This estimates the new applied bias based on linear extrapolation of the last two values. 
	//It speeds up convergenge as it avoids one loop of applied bias calculations. New applied bias is still actually determined and updates, so errors dont compound 
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
	while (std::abs(potAtReference) > 0.0001) 
	{
		//electricFieldStart is the initial guess for a boundary condition because the real boundary conditions are impossible to apply directly (assumes that the potential drops linearly over the cell at the start)
		//This boundary condition will then be updated based on the potential profile that is calculated
		//This process is repeated untal a potential profile is generated that satisfies the real boundary conditions (pot at negative electrode == applied potential && pot at reference electrode == 0)
		double electricFieldStart{ m_newAppliedBias / m_thickness };
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
			differenceEF = (m_newAppliedBias - electricFieldIntegral) / m_newAppliedBias;
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
		m_newAppliedBias += potAtReference / m_referencePositionRelative;
	}
	m_oldAppliedBias = m_appliedBias;
	m_appliedBias = m_newAppliedBias; //Update the biases for the next round
}

void DisCell::inspectPotentialODE(std::ofstream& inspectionFile)
{
	//std::cout << m_dxf << '\n' << m_dxs1 << '\n' << m_dxs2 << '\n';
	//This function takes in the reference to ePotentialArray, where the electrostatic potential profile and its 2 space derivates reside
	//It updates the array directly based on the input and thus does not need to return any values. vBias is also updated live to speed up convergence in the next steps

	inspectionFile << "-----Inspection of potential ODE started-----\n";
	inspectionFile << "Current applied bias: " << m_electrostatic[es_potential][0] << '\n';
	unsigned int totalOuterCycles{ 0 };
	unsigned int totalInnerCycles{ 0 };

	m_newAppliedBias = (2 * m_appliedBias - m_oldAppliedBias);
	//calculate the second space derivative of the electrostatic potential for every cell based on the Poisson equation
	for (array_type::size_type i{ 1 }; i < m_interfacePoint + 1; ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantFilm * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i] - m_concentrations[carrier_anions][i]);
	}

	for (array_type::size_type i{ m_interfacePoint + 1 }; i < (m_size - 1); ++i)
	{
		m_electrostatic[es_secondDerivative][i] = m_poissonConstantSolution * (-m_concentrations[carrier_electrons][i] + m_concentrations[carrier_cations][i] - m_concentrations[carrier_anions][i]);
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
	m_appliedBias = m_newAppliedBias; //Update the biases for the next round
	inspectionFile << "Outer cycles: " << totalOuterCycles << '\n';
	inspectionFile << "Inner cycles: " << totalInnerCycles << '\n';
	inspectionFile << "-----End of this inspection-----\n\n\n";
}

void DisCell::initializeConcentrations()
{
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
	/*
	Fitting the electron DOS correction factor
		This one is a doozy and it does improve the fit with experiment noticibally but not by that much
		It works as follows :
	1) We fill up the DOS with moreand more electrons.If there are more electrons in a certain part of the QD film, they will be in higher energy levels.
	2) So the Fermi level in parts with more electrons lies higher, giving rise to additional force(besides normal diffusion) that pushes the electrons towards regions with less concentration.
	3) We want to know the difference in Fermi level between each cell, so that we can easily caculate the "total" energy level of electrons(= Fermi level + electrostic potential)
	4) We then use that total energy level instead of just the electrostatic potential level in a cell for the drift - diffusion equations
	5) Hope you are still fowllowing this.
	6) So we want to know the Fermi level in each cell, which we can get by integrating the DOS, but this is computationally expensive.If we need to do this for every cell every timestep, we slow down the simulation A LOT
	7) So instead we hope to find a function that can be used to calculate the Fermi level directly from the electron concentration.
	8) Biggest assumption : DOS function is sqrt with energy.This is almost true for ZnO, not so much for CdSe maybe but it still fits rather OK.
	9) So then doing the math the integral of the DOS function should be a function of E ^ (3 / 2)
	10) So electron concentration scales with E ^ (3 / 2)
	11) And we can fit the reverse function to the electron concentration to obtain the Fermi energy
	12) Fermi level ~electron concentration ^ (2 / 3)
	13) What we really want is the difference between 2 energy levels though, so in the simulator, we use the differential function dE / dn ~n ^ (-1 / 3))
	14) So you will find in the simulator that dE = fit factor * dn / n ^ (1 / 3).Which is really complicated but it is also the fastest calculation possible.
	IF YOU FEEL THIS IS TOO COMPLICATED JUST REPLACE THE FUNCTION CALLS FOR negativeCurrente WITH THE NORMAL negativeCurrent FUNCTION!!
	*/
	return (-concentrationLeft * (electricField - m_electronEnergyFactor * (concentrationLeft - concentrationRight) / pow((concentrationLeft > 0) ? concentrationLeft : 1, 0.3333333333)) + eCon * (concentrationLeft - concentrationRight)) * curCon;
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
	for (array_type::size_type i{ 1 }; i < m_interfacePoint -3; ++i)
	{
		{
			m_currents[carrier_electrons][i] = negativeCurrente(m_concentrations[carrier_electrons][i - 1], m_concentrations[carrier_electrons][i],
				m_currentConstantElectrons, m_electrostatic[es_electricField][i - 1], m_energyConvertxf);
		}
	}

	//These lines prevent the simulation from crashing at nearly the end, when the extraction takes place. 
	//They prevent the electron concetration from going under zero near the interface during extraction.
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
		m_currents[carrier_anions][i] = negativeCurrent(m_concentrations[carrier_anions][i - 1], m_concentrations[carrier_anions][i],
			m_currentConstantAnionsFilm, m_electrostatic[es_electricField][i - 1], m_energyConvertxf);
	}

	for (array_type::size_type i{ m_interfacePoint-3 }; i < m_interfacePoint + 1; ++i)
	{
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
	m_currents[carrier_cations][m_interfacePoint + 1] = positiveCurrent(m_concentrations[carrier_cations][m_interfacePoint], m_concentrations[carrier_cations][m_interfacePoint + 1] * m_QDFillFactor,
		m_currentConstantCationsFilm2, m_electrostatic[es_electricField][m_interfacePoint], m_energyConvertxs1);
	m_currents[carrier_anions][m_interfacePoint+1] = negativeCurrent(m_concentrations[carrier_anions][m_interfacePoint], m_concentrations[carrier_anions][m_interfacePoint + 1] * m_QDFillFactor,
		m_currentConstantAnionsFilm2, m_electrostatic[es_electricField][m_interfacePoint], m_energyConvertxs1);


	//update total electrons that have entered since the last getCurrent() call:
	m_currentCumulative -= m_currents[carrier_electrons][2] + m_currents[carrier_anions][2] - m_currents[carrier_cations][2];
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
} //increases the applied bias by vbiaschange

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
