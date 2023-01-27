#include <fstream>
#include <string>
#include <algorithm>
#include <array>
#include <iostream>
#include "Config.h"

//Loads everything from the config file, and the DOS and oxidation parameters
void config(settings_array& settings, Mode& mode, ScanMode& scanmode, DOS_array& DOS, std::string configLocation)
{
    //Loading the settings and modes
    std::ifstream cFile(configLocation + "Config.txt");
    if (cFile.is_open())
    {
        std::string line;
        std::size_t index{ 0 };
        while (getline(cFile, line))
        {
            line.erase(std::remove_if(line.begin(), line.end(), isspace), //Parsing the config file
                line.end());
            if (line.empty() || line[0] == '#')
            {
                continue;
            }
            auto delimiterPos = line.find("="); //still parsing, finding the actual numbers
            if (!index)
                mode = static_cast<Mode>(std::stoi(line.substr(delimiterPos + 1))); //Selects the class/mode
            else if (index == 1)
                scanmode = static_cast<ScanMode>(std::stoi(line.substr(delimiterPos + 1)));//Selects the scanmode
            else
            {//Loads the rest of the settings as doubles
                auto value = std::stod(line.substr(delimiterPos + 1));
                settings[index - 2] = value;
            }
            index++;
        }
    }

    else
    {
        std::cerr << "Config.txt could not be opened for writing.\n";
        throw -1;      
    }

    //Loading the DOS
    std::ifstream DOSFile(configLocation + "DOS.csv");
    if (DOSFile.is_open())
    {
        DOS_array::size_type index{ 0 };
        while (DOSFile)
        {
            DOSFile >> DOS[index++];
        }
    }

    else
    {
        std::cerr << "DOS.csv could not be opened for writing.\n";
        throw - 1;
    }

    if (scanmode == ScanMode::m_CurrentTime)
    {
        settings[s_startVoltage] = settings[s_appliedBias];
    }

    settings[39] = DOS[2];      //This is the "secret" dE/dn factor that allows the correction for electron energy in the drift-diffusion equations
                                //It has to be moved from the DOS to the settings to avoid loading an extra array upon class initialization
}