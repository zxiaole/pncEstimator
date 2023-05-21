# Estimate the Particle Number Concentration Using Machine Learning

We have developed a machine learning (ML) model to estimate the particle number concentration (PNC), which is mainly induced by the ultrafine particles below 100 nm. The ML model was trained using the unique long-term observations of PNC in Switzerland by the NABEL network.

## Downscale the resolution of CAMS
* [downScaleConc.m](/src/downScaleConc.m): we prepare the data for training a machine learning model to estimate the NOx concentration at the monitoring sites of NABEL. There are hourl measurment of NOx concentrations at 16 sites in NABEL network. The data from air quality model in the CAMS database is only 10 km. We try to develop a model to estimate the measurement based on the data of the grid where the sites are located in the data from CAMS, in conjuction with  other features, e.g. wind, solar radiation and others. 

