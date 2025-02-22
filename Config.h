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
	m_discontinuous,
	m_discontinuousSideReaction,
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
	s_amountOfFilmCells,
	s_amountofInterfaceCells,
	s_interfaceResolution,
	s_dt,
	s_startVoltage,
	s_stopVoltage,
	s_scanRate,
	s_voltageIncrement,
	s_appliedBias,
	s_runTime,
	s_electronMobility,
	s_cationMobilityFilm,
	s_cationMobilitySolution,
	s_anionMobilityFilm,
	s_anionMobilitySolution,
	s_redoxSpeciesConcentration,
	s_ionConcentration,
	s_QDspacefill,
	s_temperature,
	s_epsilonrFilm,
	s_epsilonrSolution,
	s_LUMO,
	s_negativeElectrodeWF,
	s_kineticsType,
	s_contaminantConcentration,
	s_contaminantMobilityFilm,
	s_contaminantMobilitySolution,
	s_BV_E0,
	s_BV_k,
	s_kIrriversible,
	max_settings
};

namespace physics //contains all natural constants, access through physics::
{
	constexpr double k{ 1.380649e-23 };
	constexpr double eps0{ 8.854188e-12 };
	constexpr double h{ 6.62607004e-34 };
	constexpr double q{ 1.60217662e-19 };
	constexpr double F{ 96485 };
	constexpr double R{ 8.31446 };
}

void config(settings_array& settings, Mode& mode, ScanMode& scanmode, DOS_array& DOS,std::string configLocation);
