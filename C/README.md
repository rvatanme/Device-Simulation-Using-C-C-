# Numerical solution of DD current equations using C Programming language
Here the code written in C has been provided to simulate the IV curve of a one-dimensional (1d) abrupt Si PN diode based on DD model.

The code has been devided to several parts and the comments are added to the code where ever the explanation is needed.

The code takes 2.0 s to complete the job. During running, several text files are produced that contains related physical quantities of the simulated 1d PN diode that can be plotted by any desire plotting software. 

The generated files are as follows:

1) "els_plot.txt"     --->     Contains the Potential (V), Electric Field (V/cm), Electron Densities (1/cm^3), Hole Densities (1/cm^3), Total Charge Density (C/cm^3),	Conduction Band Energy (eV) vs x (um) at the equilibrium when no external potential is applied on the diode. The graphs plotted by libreoffice have been provided in the current folder. 
2) "quasi.txt" ---> contains the conduction, valence and interinsic band energy (eV) and quasi Fermi level (eV) vs x (um) at the applied voltage of 0.625 V.
3) "pote.txt"  ---> contains potential profile (V) vs position (um) at the applied voltage of 0.625 V.
4) "elec.txt"  ---> contains electric field profile (MV/cm) vs position (um) at the applied voltage of 0.625 V.
5) "car_dens.txt"  ---> contains electron and hole densities (/cm3) vs position (um) at the applied voltage of 0.625 V.
6) "curr_tot_ele_hol.txt"  ---> contains Total Current Density (Amp/cm^2)	Electron Current Density (Amp/cm^2)	Hole Current Density (Amp/cm^2) vs Applied Voltage (V) at the applied voltage of 0.625 V.
7) "curr_pos.txt"  ---> contains Total Current Density (Amp/cm^2) vs position (um) at the applied voltage of 0.625 V.
