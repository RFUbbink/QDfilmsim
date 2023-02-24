# QDfilmsim by Reinout Ubbink
Drift-diffusion simulator of electrochemistry on QD films

Welcome.
You will find both the source code and a windows executable. 
Source code contains its own comments.
Executable needs a config file in the working directory (Config.txt; example file included), in which simulation parameters are set. 
Altering anything else than the numbers after the '=' sign may cause crashes, as the parsing algorithm is very primitive. 
Executable also needs a density of states function file (DOS.csv; example file included) in the working directory. 
Python script DOSdrawer.py can be used to produce a DOS file in the correct formatting. Alter the functions to adjust the DOS output.
Executable will output 'outputcpp.csv', containing the CV or time/current data and 'Endfile.csv, containing the entire final state of the simulation.
Executable will also output 'Midfile#.csv' after every 0.1 V has passed in a CV. This file contains the entire state of the simulation, including currents, potential and concentrations.
