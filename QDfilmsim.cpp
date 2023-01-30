/*
Welcome to the QD film electrochemistry simulator source code.
Good luck in traversing it! If you need help reading/learning coding concepts of C++, go to: learncpp.com. All the classes are explained in their corresponding HEADER files. 
Some classes are legacy content, but I did not feel like taking them out in case they contained some usefull bits of code that I would need later.
I will take you through the main function in comments. 
*/

#include <iostream>
#include <math.h>
#include <vector>
#include <array>
#include <cstddef>
#include <chrono>
#include <fstream>
#include <string>
#include "Electrochemistry.h"
#include "Xcontamination.h"
#include "RIreaction.h"
#include "Discontinuous.h"
#include "DisXcon.h"
#include "OxidationGM.h"
#include "DisOxGM.h"
#include "NoFilm.h"
#include "Config.h"


void currentToFile(std::ofstream& outf, std::vector<double>& voltage, std::vector<double>& electronCurrent, std::vector<double>& leakCurrent)
{
	//This fucntion saves the voltage and current output to a file as a csv. That is basically the CV as you would get it from a potentiostat.
	for (array_type::size_type i{ 0 }; i < electronCurrent.size(); ++i)
	{
		outf << voltage[i] << '\t' << electronCurrent[i] << '\t' << leakCurrent[i] << '\n';
	}
}

settings_array settings; //settings available in global space for convenience.
DOS_array DOS;

template <typename T> //Allows the multiple classes (Cell, Discell, etc. to use the same basic running functions such as runIV)
void continueIV(T& cell, int& tries, std::string saveDirectory)
{
	//Kind of a legacy function that can be used to run another cv after finishing the first one. It is prone to crashing though, as it has not been tested for a lot of different initial states.
	//See runIV for more details
	const int cyclesPerIncr{ static_cast<int>(settings[s_voltageIncrement]
						/ settings[s_scanRate] / settings[s_dt]) };	//the number of cycles after which the voltage is incremented
	const int numberOfSteps{ 2 * static_cast<int>(std::abs(settings[s_stopVoltage]
										- settings[s_startVoltage]) / settings[s_voltageIncrement]) };	//the number of steps that are taken (factor 2 for forward/backward)
	std::vector<double> voltage(numberOfSteps);
	std::vector<double> electronCurrent(numberOfSteps);					//array to store the current over time in			
	std::vector<double> leakCurrent(numberOfSteps);					//array to store the current over time in			

	int recordTime{ 45 };
	if (settings[s_recordingVoltage])
		recordTime = static_cast<int>((-settings[s_recordingVoltage] + settings[s_startVoltage]) / settings[s_voltageIncrement]);


	for (int V{ 0 }; V < numberOfSteps; ++V)
	{
		for (int t{ 0 }; t < cyclesPerIncr; ++t)
		{
			cell.calculatePotentialProfile();
			if (!(t % 10))
				cell.injectElectrons(DOS);
			cell.calculateCurrents();
			cell.updateConcentrations();
		}

		//record current density at the end of the voltage step (converted to A/cm2) 
		voltage[V] = cell.getVoltage();
		electronCurrent[V] = cell.getCurrent()/cyclesPerIncr;
		leakCurrent[V] = cell.getLeakCurrent();

		std::cout << V << '\n'; //report to console the amount of steps taken
		//manage the applied bias by incrementing
		if (V < numberOfSteps / 2)
			--cell;

		else
			++cell;

		if (V == recordTime)
		{
			std::ofstream midf{ saveDirectory + "\\Midfile.csv", std::ios::trunc };
			if (!midf)
			{
				std::cerr << "Midfile.csv could not be opened for writing.\n";
				continue;
			}
			else
				cell.midSave(midf);
		}
	}


	std::ofstream outf{ saveDirectory + "\\IV_continued_" + std::to_string(++tries) + ".csv", std::ios::trunc };

	if (!outf)
	{
		std::cerr << "Something went wrong in writing the results.\n";
	}
	else
	{
		currentToFile(outf, voltage, electronCurrent, leakCurrent);
	}
}
template <typename T>
void runIV(T& cell, std::string saveDirectory)
{
	//The main function that is important, takes the Cell object and givs it the necessary orders and does bookkeeping.
	//For information on the memberfunctions of cell, see the corresponding object files
	const int cyclesPerIncr{ static_cast<int>(settings[s_voltageIncrement]
						/ settings[s_scanRate] / settings[s_dt]) };	//the number of cycles after which the voltage is incremented
	const int numberOfSteps{ 2 * static_cast<int>(std::abs(settings[s_stopVoltage]
										- settings[s_startVoltage]) / settings[s_voltageIncrement])};	//the number of steps that are taken (factor 2 for forward+backward)
	
	std::ofstream outf{ saveDirectory + "\\Outputcpp.csv", std::ios::trunc }; //set up some files for saving
	std::string title{ saveDirectory + "\\Midfile0" };
	int recordCounter{ 1 };

	cell.initializeConcentrations(settings[s_sontaminantConcentration]); 

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now(); //to mesure performance
	//Here starts the main cycle. It will do voltage steps, and within them a number of time steps.
	for (int V{ 0 }; V < numberOfSteps; ++V)
	{

		for (int t{ 0 }; t < cyclesPerIncr; ++t)
		{
			//The 4 steps that are taken each cycle. Function speak for themselves I guess. 
			cell.calculatePotentialProfile();
			if (!(t % 10))
			{
				cell.resetInjection(); //Need this one to avoid weird stuff at the interface (negative concentrations etc.) it doesnt affect the results.
				cell.injectElectrons(DOS); //Injection is limited to once every 10 steps for performance, since an integral over the DOS needs to be performed every time.
				//Since the timesteps are short, this is not an issue, the concentration of electrons is nearly constant over this interval
			}
			cell.calculateCurrents();
			cell.updateConcentrations();
		}
		
		double current{ cell.getCurrent() }; //record current density at the end of the voltage step (converted to A/cm2) 
		if (isnan(current))
		{
			//Some very primitive crash detection. The simulation can still crash silently when the ODEsolver for the potentail profile fails to converge
			std::cout << "Simulation crash (current isNaN)\n";
			std::ofstream crashf{ saveDirectory + "\\crashfile.csv", std::ios::trunc };
			cell.midSave(crashf);
		}
		outf << settings[s_startVoltage] - settings[s_voltageIncrement]* ((V < (numberOfSteps / 2)) ? V : (numberOfSteps - V)) << '\t' << current / cyclesPerIncr << '\t' << cell.getLeakCurrent()  << '\t' << cell.getVoltage() << '\n';

		std::cout << V << '\n'; //report to console the amount of steps taken

		//midsim recording of the state of the cell, reports concentrations, currents, and electrostatic data.
		if (V%100 == 90)
		{

			std::ofstream midf{ title + std::to_string(recordCounter++) + ".csv", std::ios::trunc };
			if (recordCounter > 9)
				title = saveDirectory + "\\Midfile";
			if (!midf)
			{
				std::cerr << "Midfile.csv could not be opened for writing.\n";
				continue;
			}
			else
				cell.midSave(midf);
		}

		//manage the applied bias by incrementing
		if (V < (numberOfSteps / 2))
			--cell;

		else
			++cell;
	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "\nTime elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << '\n';

	bool continueScan{false};

	while (continueScan)
	{//To ask the user for another scan
		int tries{ 0 };
		std::cout << "Another? [y/n]: ";
		char input{};
		std::cin >> input;
		if (std::cin.fail())
		{
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		else if (input == 'y')
		{
			continueIV(cell, tries, saveDirectory);
		}
		else if (input == 'n')
		{
			continueScan = false;
		}
	}
}
template <typename T>
void runCurrentTime(T& cell, std::string saveDirectory)
{//Very similar to the runIV fucntion except that is does not increment the voltage but instead applies a constant votlage for the durration of the simulation. 
	long long int numberOfSteps{ static_cast<long long int>(settings[s_runTime] / settings[s_dt]) };	//the number of steps that are taken (factor 2 for forward/backward)
	while (numberOfSteps % 100)
		++numberOfSteps;
	const long long int resolution{ 10000 };
	std::vector<double> time(resolution);
	std::vector<double> electronCurrent(resolution);					//array to store the current over time in (no need for resolution beyond 10000)	
	std::vector<double> leakCurrent(resolution);					//array to store the current over time in (no need for resolution beyond 10000)		

	const int divisor{ static_cast<int>(numberOfSteps / resolution) };

	cell.initializeConcentrations(settings[s_sontaminantConcentration]);

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	for (long long int t{ 0 }; t < numberOfSteps; ++t)
	{
		cell.calculatePotentialProfile();
		if (!(t % 10))
		{
			cell.resetInjection();
			cell.injectElectrons(DOS);
		}
		cell.calculateCurrents();
		cell.updateConcentrations();
		if (!(t % divisor))
		{
			time[static_cast<int>(t / divisor)] = t * settings[s_dt];
			electronCurrent[static_cast<int>(t / divisor)] = cell.getCurrent(); //record current density at the end of the voltage step (converted to A/cm2) 
			leakCurrent[static_cast<int>(t / divisor)] = cell.getLeakCurrent(); //record current density at the end of the voltage step (converted to A/cm2) 
			if (!((100 * t) % numberOfSteps))
				std::cout << 100 * t / numberOfSteps << " % complete.\n";
		}
	}
	
	std::ofstream endf{ saveDirectory + "\\Endfile.csv", std::ios::trunc };
	if (!endf)
	{
		std::cerr << "Endfile.csv could not be opened for writing.\n";
	}
	else
		cell.midSave(endf);

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "\nTime elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << '\n';

	//time to write the results to a file for Phyton analysis
	std::ofstream outf{ saveDirectory + "\\Outputcpp.csv", std::ios::trunc };

	if (!outf)
	{
		std::cerr << "Outputcpp.csv could not be opened for writing.\n";
	}
	else
		currentToFile(outf, time, electronCurrent, leakCurrent);
}

int main(int argc, char* argv[])
{
	//Main simply loads the config data, then selects a class and runs a function based on the settings. 
	//Note: the sim requires a simple file structure, where the executable is in the same folder as the config file. 
	Mode mode{ Mode::m_standard };
	ScanMode scanmode{ ScanMode::m_IVcurve };
	std::string configLocation(argv[0]); //First argument is always the location of the program itself
	configLocation = configLocation.substr(0, configLocation.size() - 13); //We cut off the executable part (13 characters, 15 for Linux) to find the config file in the SAME DIRECTORY AS THE EXECUTABLE
	try
	{
		config(settings, mode, scanmode, DOS, configLocation); //load settings from config file (Config.txt)
	}		
	catch (int)
	{
		std::cerr << "Failed to load config file!\n";
		return 1;
	}
	
	std::string saveDirectory(configLocation); //If no save directory is provided, default to saving at the config location. 
	if (argv[1])
		saveDirectory = std::string(argv[1]); //If save directory is provided, select it. 
	std::cout << saveDirectory << '\n';
	if (scanmode == ScanMode::m_IVcurve) //Select the operating mode
	{
		switch (mode) //Select the right object type
		{
		case Mode::m_standard:
		{
			Cell cell{ settings }; //Ready to go, lets start the loop!
			runIV<Cell>(cell, saveDirectory);
			break;
		}
		case Mode::m_xcontamination:
		{
			Xcontamination xcontamination{ settings };
			Cell& cell{ xcontamination }; //Ready to go, lets start the loop!
			runIV<Cell>(cell, saveDirectory);
			break;
		}
		case Mode::m_rireaction:
		{
			RIreaction rireaction{ settings };
			Cell& cell{ rireaction }; //Ready to go, lets start the loop!
			runIV<Cell>(cell, saveDirectory);
			break;
		}
		case Mode::m_discontinuous:
		{
			DisCell discontinuous{ settings };//Ready to go, lets start the loop!
			runIV<DisCell>(discontinuous, saveDirectory);
			break;
		}
		case Mode::m_oxidationGM:
		{
			DOS_array RR{};
			std::ifstream RRFile(configLocation + "ReactionRates.csv");
			if (RRFile.is_open())
			{
				DOS_array::size_type index{ 0 };
				while (RRFile)
				{
					RRFile >> RR[index++];	
				}
			}

			else
			{
				std::cerr << "ReactionRates.csv could not be opened for writing.\n";
				return 1;
			}
			OxidationGM oxidationGM{ settings, RR };//Ready to go, lets start the loop!
			Cell& cell{ oxidationGM };
			runIV<Cell>(cell, saveDirectory);
			break;
		}
		case Mode::m_DisOxGM:
		{
			DOS_array RR{};
			std::ifstream RRFile(configLocation + "ReactionRates.csv");
			if (RRFile.is_open())
			{
				DOS_array::size_type index{ 0 };
				while (RRFile)
				{
					RRFile >> RR[index++];
				}
			}

			else
			{
				std::cerr << "ReactionRates.csv could not be opened for writing.\n";
				return 1;
			}
			DisOxGM disOxGM{ settings, RR};//Ready to go, lets start the loop!
			DisCell& discell{ disOxGM };
			runIV<DisCell>(discell, saveDirectory);
			break;
		}
		case Mode::m_NoQDFilm:
		{
			NoFilmCell noFilmCell{ settings};//Ready to go, lets start the loop!
			runIV<NoFilmCell>(noFilmCell, std::string(argv[0]));
			break;
		}
		}
	}

	else if (scanmode == ScanMode::m_CurrentTime)
	{
		switch (mode)
		{
		case Mode::m_standard:
		{
			settings[s_startVoltage] = settings[s_appliedBias];
			Cell cell{ settings };
			runCurrentTime<Cell>(cell, saveDirectory);
			break;
		}
		case Mode::m_xcontamination:
		{
			settings[s_startVoltage] = settings[s_appliedBias];
			Xcontamination xcontamination{ settings };
			Cell& cell{ xcontamination }; //Ready to go, lets start the loop!
			runCurrentTime<Cell>(cell, saveDirectory);
			break;
		}
		case Mode::m_rireaction:
		{
			RIreaction rireaction{ settings };
			Cell& cell{ rireaction }; //Ready to go, lets start the loop!
			runCurrentTime<Cell>(cell, saveDirectory);
			break;
		}
		case Mode::m_discontinuous:
		{
			DisCell discontinuous{ settings };//Ready to go, lets start the loop!
			runCurrentTime<DisCell>(discontinuous, saveDirectory);
			break;
		}
		default:
		{
			std::cerr << "Mode not yet implemented for current-time experiments!\n";
			return 1;
		}
		}

	}
	argc; //Throw argc away to avoid warnings
	return 0;
}

