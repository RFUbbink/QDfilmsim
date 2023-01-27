#pragma once

using array_type = std::vector<std::vector<double>>;

namespace config // contains all the configuration settings that should be moved to a config file
{
	//spacetime parameters
	inline constexpr double filmThickness{ 100e-9 };		//QD film thickness in m
	inline constexpr double referencePosition{ 200e-9 };	//position of the reference electrode relative to the working electrode in m
	inline constexpr double cellThickness{ 300e-9 };		//distance between the working- and counterelectrode in m
	inline constexpr int	amountOfCells{ 150 };			//amount of cells in the EC cell
	inline constexpr double workingArea{ 1e-4 };			//area of the working electrode in m2
	inline constexpr double counterArea{ 1e-4 };			//area of the counterelectrode in m2
	inline constexpr double dt{ 1e-10 };					//time step size is s

	//scan characteristics
	inline const double startVoltage{ -0.1 };				//starting and ending voltage of the scan
	inline const double stopVoltage{ -0.6 };					//turning point voltage of the scan
	inline const double scanRate{ 10000 };					//scan rate in Vs-1
	inline const double increment{ 0.01 };					//Voltage increment in V

	//mobilities of the carriers in m2V-1s-1
	inline constexpr double mobilityElectrons{ 5e-8 };
	inline constexpr double mobilityCationsFilm{ 5e-8 };
	inline constexpr double mobilityCationsSolution{ 5e-8 };
	inline constexpr double mobilityAnionsFilm{ 5e-8 };
	inline constexpr double mobilityAnionsSolution{ 5e-8 };

	//misc
	inline constexpr double denstityOfStates{ 3e25 };		//effective electron density of states in m-3
	inline constexpr double saltConcentration{ 1.25e25 };	//initial ion concentration in solution in m-3
	inline constexpr double QDdensity{ 0.5 };				//the fraction of space filled by the QD, the rest is solvent. Used to calculate the ion concentration in the film
	inline constexpr double temperature{ 300.0 };			//in K
	inline constexpr double epsilonrFilm{ 4.0 };			//relative dielectric constant in a quantum dot film
	inline constexpr double epsilonrSolution{ 37.0 };		//relative dielectric constant in an ACN solution

	//band characteritics
	inline constexpr double LUMO{ 0.5 };					//energy level of the conduction band relative to vacuum
	inline constexpr double negativeElectrodeWorkFunction{ 0 }; //work function of the negative electrode compared to vacuum
}

namespace phys //contains all natural constants, access through phys::
{
	inline constexpr double k{ 1.380649e-23 };
	inline constexpr double eps0{ 8.854188e-12 };
	inline constexpr double h{ 6.62607004e-34 };
	inline constexpr double q{ 1.60217662e-19 };
}


