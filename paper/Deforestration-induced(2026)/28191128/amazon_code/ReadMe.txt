Folder Content:
- Historic climate change data (historic_climate_change_data)
- Deforestation/Land-Use-Change Data for SSP scenarios. For deforestation data by Soares-Filho et al (2013), visit: https://doi.org/10.3334/ORNLDAAC/1153 (deforestation_data)
- average input files for MAP, MCWD and the moisture recycling network (average_network)
- Main Simulation files (main_simulations)
- ensemble construction (probabilistic_ensemble)
- software to solve nonlinear interacting differential equations (pycascades)


How to start ensemble of simulations (in the folder main_simulations):
1) Starte the file "drought_unstable_amazon.py" to iniate tipping experiments using:
	(i) different initial conditions on the adaptive capacity of the Amazon forest using the ensemble stored in the folder "probabilistic_ensemble"
	(ii) Different climate change scenarios using the scenarios stored in the folder "average_network"
	(iii) If you also aim to include deforestation scenarios, download the respective deforestation dataset and insert the deforested parts of the forest in the variable >> init_state.fill(-1) instead of "-1"
2) Collect the results and map out tipping likelihoods for the different scenarios 
3) Add a tipping reason analysis with the file "drought_unstable_amazon_tipping_reason.py"



Remarks:
1) Not all original data is stored in these folders due to space limitations but their original source is given in the manuscript's data (and code) availability statement.
2) In case of questions regarding the code or requests for data, please be invited to contact the corresponding authors of this study. We are happy to help and get you ready to simulate your own scenarios.
