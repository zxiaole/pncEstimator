# Estimate the Particle Number Concentration Using Machine Learning

We have developed a machine learning (ML) model to estimate the particle number concentration (PNC), which is mainly induced by the ultrafine particles below 100 nm. The ML model was trained using the unique long-term observations of PNC in Switzerland by the NABEL network.

## Downscale the NOx concentrations with low-resolution of CAMS
* [downScaleConc.m](/src/preProcessing/downScaleConc.m): we prepare the data for training a machine learning model to estimate the NOx concentration at the monitoring sites of NABEL. There are hourl measurment of NOx concentrations at 16 sites in NABEL network. The data from air quality model in the [CAMS database](https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-europe-air-quality-forecasts?tab=form) is only 10 km. We try to develop a model to estimate the measurement based on the data of the grid where the sites are located in the data from CAMS, in conjuction with  other features, e.g. wind, solar radiation and others. The meteo data are from the [ERA5-Land hourly data](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=overview). The traffic data are from the [Open Transport Map](https://hub.plan4all.eu/otm), which provides the [database](https://opentransportmap.info/download) for the European countries.

* [downScaleModel.ipynb](/src/processing/downScaleModel.ipynb): The data generated above will be imported into the python code to train a machine learning model. The code can be used to (1) generate the ML model using the NABEL data; (2) generate the map with hight resolution for the whole Switzerland using the trained model and the gridded data.
* [swissPNCDistribution_downScale.m](/src/postProcessing/swissPNCDistribution_downScale.m): the code is used to generate the downscaled map as shown below
 <img src="/img/DownScale_NOX.png">
