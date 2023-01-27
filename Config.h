/*
Mainly contains some enums, to do with loading the config data. Made the enums to make it easier and not get confused when adding new modes/parameters
Also contains the physical constants (phys)
Also contains the configuration function itself, which looks at the config file and basically puts it all in the RAM.
*/

#pragma once
#include <vector>
#include <array>
#include <string>

using settings_array = std::array<double, 40>;
using DOS_array = std::array<double, 402>;
//The RR array's first value is the maximum electron concentration. The array will then assume a linear electron concentration
using array_type = std::vector<std::vector<double>>;

enum class Mode //For the program to know which mode it should take. Each mode has a corresponding class
{
	m_standard,
	m_xcontamination,
	m_rireaction,
	m_discontinuous,
	m_oxidationGM,
	m_DisOxGM,
	m_NoQDFilm,
};

enum class ScanMode //The scanmode is selected in the main() function
{
	m_IVcurve,				//Normal current versus voltage curve
	m_CurrentTime,			//Current vs time 
};

enum settings //enum that is used for sorting the data from the config file/array
{
	s_filmThickness,
	s_refPosition,
	s_cellThickness,
	s_amountOfCells,
	s_workingElectrodeArea,
	s_counterElectrodeArea,
	s_dt,
	s_startVoltage,
	s_stopVoltage,
	s_recordingVoltage,
	s_scanRate,
	s_voltageIncrement,
	s_appliedBias,
	s_runTime,
	s_electronMobility,
	s_cationMobilityFilm,
	s_cationMobilitySolution,
	s_anionMobilityFilm,
	s_anionMobilitySolution,
	s_densityOfStates,
	s_ionConcentration,
	s_QDspacefill,
	s_temperature,
	s_epsilonrFilm,
	s_epsilonrSolution,
	s_LUMO,
	s_negativeElectrodeWF,
	s_contaminantMobility,
	s_sontaminantConcentration,
	s_maximumElectronConcentration,
	s_maximumIonConcentration,
	s_k1,
	s_E1,
	s_k2,
	s_E2,
	s_maxTrapConcentration,
	max_settings
};

namespace precompiled
{
	inline constexpr int	amountOfCells{ 490 };			//amount of cells in the EC cell, this is precompiled because of legacy, but could easily not be now that I use classes etc. 
	inline constexpr int	InterfaceCells{ 60 };			//amount of cells that define the interface region
}

namespace phys //contains all natural constants, access through phys::
{
	inline constexpr double k{ 1.380649e-23 };
	inline constexpr double eps0{ 8.854188e-12 };
	inline constexpr double h{ 6.62607004e-34 };
	inline constexpr double q{ 1.60217662e-19 };
	inline constexpr double F{ 96485 };
	inline constexpr double R{ 8.31446 };
}

void config(settings_array& settings, Mode& mode, ScanMode& scanmode, DOS_array& DOS,std::string configLocation);
