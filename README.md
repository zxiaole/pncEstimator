# Estimate the Particle Number Concentration Using Machine Learning

We have developed a machine learning (ML) model to estimate the particle number concentration (PNC), which is mainly induced by the ultrafine particles below 100 nm. The ML model was trained using the unique long-term observations of PNC in Switzerland by the NABEL network.

## Downscale the NOx concentrations with low-resolution of CAMS
* [downScaleConc.m](/src/preProcessing/downScaleConc.m): we prepare the data for training a machine learning model to estimate the NOx concentration at the monitoring sites of NABEL. There are hourl measurment of NOx concentrations at 16 sites in NABEL network. The data from air quality model in the [CAMS database](https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-europe-air-quality-forecasts?tab=form) is only 10 km. We try to develop a model to estimate the measurement based on the data of the grid where the sites are located in the data from CAMS, in conjuction with  other features, e.g. wind, solar radiation and others. The meteo data are from the [ERA5-Land hourly data](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=overview)

