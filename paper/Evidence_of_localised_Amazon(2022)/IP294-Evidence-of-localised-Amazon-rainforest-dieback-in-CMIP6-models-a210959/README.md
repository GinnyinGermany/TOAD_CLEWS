# Evidence-of-localised-Amazon-rainforest-dieback-in-CMIP6-models
## Project Description
The project aimed to design an algorithm which would detect abrupt shifts in CMIP6 data and use this to observe potential tipping points in vegetation carbon in the Amazon rainforest within different CMIP6 models. Further analysis looked into a potential indicator of risk for grid points in the Amazon which were detectd as having an abrupt shift

The scripts in this repository take global cveg and tas data from CMIP6 models and interplolates it onto a 1 degree world grid, then interpolates for a specific region (Amazon). An algorithm is then used on the interpolated data, along with control run data, to detect abrupt shifts analogous to tipping events in the vegetation carbon in the Amazon rainforest. The detected abrupt shifts may then be analysed to observe how they evolve with time and CO2 levels, how many abrupt shifts are observed in the region and to find a metric for the risk associated with an abrupt shift occuring for these grid points. 

## Script Run Order
1. World_interpolation / World_interpolation_tas
2. Regional_interpolation / Regional_interpolation_tas
3. windows_gradient_picontrol
4. abrupt_shift_detection_statistics
5. monthly_statistical_analysis
6. All other scripts may now be run

## Authors and Acknowledgment 
Isobel Parry, Paul Ritchie and Peter Cox
