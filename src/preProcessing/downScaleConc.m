%% Aug. 11, 2021, Xiaole Zhang
% The original resolution of CAMS air quality data is 10 km
% We try to downscale the data to 1 km resolution using machine learning
% method with NABLE measurements and features.
% Downscale the concentration from CAMS based on the traffic volume, population density and other features
% from https://ghsl.jrc.ec.europa.eu/ghs_pop2019.php

%% Sep. 8, 2021, Xiaole Zhang
% Use the air pollution data from CAMS as indicator;
% Develop a model to estimate the deviation of the measurements from CAMS;
% Parameters (refer to the section "%% combine all the data" below):
% (1) low-resoultuion (10km) grid data from CAMS
% (2) solar radiation (W/m*2)
% (3) surface temperature at 2m (celsius degree)
% (4) precipitation (mm/h)
% (5) relative humidity at 2m (%, estimated by dew temperature)
% (6) wind speed at 10m (m/s)
% (7) road traffic volume
% (8) pop density;
% (9) hour;
% (10) month;
% (11) weekday

% the target variable is the ratio between the observations from NABLE
% statoins and the low-resolution CAMS data at the corresponding location

% we use the data to train a ML model (e.g. random forest) to estimate the
% target variable based on the 11 parameters metioned above.

% this code prepares all the data needed for training the ML model.
% later the data will be feeded into the python code to train the model.

% predictionFlag = 1 save the data for the downscale estimation for the
% whole Switzerland;
% predictionFlag = 0 only save the data with NABEL measurements at the
% sites for 2019 and 2020

%% May 19, 2023, Xiaole Zhang
% Organize the codes for github

%% clean the workspace
clear

%% config
dataFolder = '../data';

% resolution for the downscaling
Delta_x=0.01; % in degree, about 1 km
Delta_y=0.01;

% We only donscale NOX or NO2.
% The high-resolution NOX and NO2 can be used as an indicator to estimate
% high resolution distribution of PNC. However, it might be more reasonble
% to also downscale other pollutants, e.g. pm2.5 and pm10. Extra work will
% be needed but the methods for downscaling would remain similar.
pollutantType = 'NOX'; % 'NOX' or 'NO2'

% prepare the data for training '0', or for prediction '1'
predictionFlag = 0;

% Years with data
% more data can be available from NABEL, but we only use the data from
% recent years, becasue the standards for emissions keep changing over the
% years. THe data from recent years can be more representative for the
% current evaluation or prediction
years = 2019:2020;

%% load the coordinates of the stations in NABEL network
% There are much more sites for NO2 or NOx than those for PNC.
nabelSiteFile = fullfile(dataFolder, 'nabelSites.mat');
load(nabelSiteFile)

siteLat = [sites{:,3}];
siteLon = [sites{:,4}];
inputFormat = 'dd.MM.yyyy HH:mm';

%% read pollutants
folder = fullfile(dataFolder,'pollutants');
pollutant = pollutantType;

for yearId = years
    filename = fullfile(folder, [pollutant num2str(yearId) '.csv']);
    
    if(yearId == years(1))
        pollData = readtable(filename,'HeaderLines',5);
    else
        pollDataTmp = readtable(filename,'HeaderLines',5);
        pollData = vertcat(pollData,pollDataTmp);
    end
end

pollData.Date_time = datetime(pollData.Date_time, 'InputFormat',inputFormat, 'TimeZone', 'Europe/Berlin');
timeMeas = hours(pollData.Date_time - pollData.Date_time(1));

oldvariables = pollData.Properties.VariableNames;
newvariables = sites(:,2);

[~,LOCB] = ismember(newvariables,oldvariables);
newTable = pollData(:,LOCB);

%% read conc
timeDate = [];
camsNOX = [];
pollutantFolder = 'CAMS';

if(~predictionFlag)
    for year = years
        firstDay = datetime(['31.12.' num2str(year-1) ' 00:00'],'InputFormat',inputFormat, 'TimeZone', 'UTC');
        formatOut = 'yyyymmdd';
        
        pollutantFile = fullfile(dataFolder, pollutantFolder, ['SingleLevel_' num2str(year-1) '1231.nc']);%SingleLevel_20190331
        lon = ncread(pollutantFile, 'longitude');
        lat = ncread(pollutantFile, 'latitude');
        [x,y] = meshgrid(lon,lat);
        
        
        for days = 0:30:365
            currentTime = firstDay+days;
            dateStr = datestr(currentTime, formatOut);
            pollutantFile = ['SingleLevel_' dateStr '.nc'];%SingleLevel_20190331
            filename = fullfile(dataFolder, pollutantFolder, pollutantFile);
            
            noxSwiss = ncread(filename, 'no2_conc')+ ncread(filename, 'no_conc');
            
            lastId = size(camsNOX,1);
            for hourId = 1:size(noxSwiss,4)
                tmp = squeeze(noxSwiss(:,:,1,hourId))';
                camsNOX(lastId+hourId, :) = interp2(x, y, tmp, siteLon, siteLat);
            end
            timeDate = [timeDate; double(ncread(filename, 'time'))/24+currentTime];
        end
    end
    % interp to the measurement time
    timeCAMS = hours(timeDate - pollData.Date_time(1));
    [vv, idV] = unique(timeCAMS);
    
    camsNOX = interp1(timeCAMS(idV), camsNOX(idV,:), timeMeas);
else
    firstDay = datetime(['31.12.' num2str(year-1) ' 00:00'],'InputFormat',inputFormat, 'TimeZone', 'UTC');
    formatOut = 'yyyymmdd';
    
    pollutantFile = ['SwissGridData/CAMS/SingleLevel_' num2str(year-1) '1231.nc'];%SingleLevel_20190331
    lon = ncread(pollutantFile, 'longitude');
    lat = ncread(pollutantFile, 'latitude');
    [x,y] = meshgrid(lon,lat);
    [xCAMS,yCAMS] = meshgrid(lon,lat);
    
    camsAllData = [];
    
    for days = 0:30:365
        currentTime = firstDay+days;
        dateStr = datestr(currentTime, formatOut);
        pollutantFile = ['SingleLevel_' dateStr '.nc'];%SingleLevel_20190331
        filename = fullfile(dataFolder, pollutantFolder, pollutantFile);
        
        if(strcmp(pollutantType, 'NOx'))
            no2Swiss = ncread(filename, 'no2_conc')+ ncread(filename, 'no_conc');
        elseif(strcmp(pollutantType, 'NO2'))
            no2Swiss = ncread(filename, 'no2_conc');
        end
        
        camsAllData = cat(3, camsAllData, squeeze(no2Swiss));
        timeDate = [timeDate; double(ncread(filename, 'time'))/24+currentTime];
    end
    
    
    timeCAMS = hours(timeDate - pollData.Date_time(1));
    
    lonNew = lon(1):Delta_x:lon(end);
    latNew = lat(1):-Delta_y:lat(end);
    [xnew, ynew] = meshgrid(lonNew, latNew);
    
    %     camsNO2 = interp2(x, y, squeeze(no2Swiss(:,:,1,1))', xnew, ynew);
end

%% read meteo
meteoRadiation = [];
meteoTemperature = [];
meteoPrecipitation = [];
meteoHumidity = [];
meteoSpeed = [];
timeMeteo = [];

if(~predictionFlag)
    for year = years
        meteoFolder = 'meteo';
        meteoFile = [num2str(year) '.nc'];%SingleLevel_20190331
        
        filename = fullfile(dataFolder, meteoFolder, meteoFile);
        ncdisp(filename);
        lon = ncread(filename, 'longitude');
        lat = ncread(filename, 'latitude');
        [x,y] = meshgrid(lon,lat);
        
        inputFormat = 'dd.MM.yyyy HH:mm';
        tmp = datetime('01.01.1900 00:00','InputFormat',inputFormat, 'TimeZone', 'UTC');
        formatOut = 'yyyymmddHH';
        currentTime = ncread(filename, 'time');
        meteoTime = tmp+double(currentTime)/24;
        
        timeMeteoTmp = hours(meteoTime - pollData.Date_time(1));
        timeMeteo = [timeMeteo; timeMeteoTmp];
        
        % 1st column: radiation
        radiation =  ncread(filename, 'ssrd');
        radiation = max(diff(radiation,1,3), 0)/3600;
        radiation = cat(3, zeros(size(radiation,1),size(radiation,2)), radiation);
        
        
        % 2ed column: temperature
        temperature = ncread(filename, 't2m')-273.15; % to celsius degree
        
        % 3rd column: precipitation
        precipitation = ncread(filename, 'tp')*1000;
        precipitation = max(diff(precipitation,1,3), 0);
        precipitation = cat(3, zeros(size(precipitation,1),size(precipitation,2)), precipitation);
        
        % 4th column: relative humidity
        dew = ncread(filename, 'd2m')-273.15; % to celsius degree
        humidity = 100*exp(17.625*dew./(243.04+dew))./exp((17.625*temperature)./(243.04+temperature));
        
        % 5th column: wind speed
        u = ncread(filename, 'u10');
        v = ncread(filename, 'v10');
        speed = sqrt(u.^2+v.^2);
        
        meteoRadiationTmp = zeros(length(timeMeteoTmp), size(sites,1));
        meteoTemperatureTmp = zeros(length(timeMeteoTmp), size(sites,1));
        meteoPrecipitationTmp = zeros(length(timeMeteoTmp), size(sites,1));
        meteoHumidityTmp = zeros(length(timeMeteoTmp), size(sites,1));
        meteoSpeedTmp = zeros(length(timeMeteoTmp), size(sites,1));
        
        for hourId = 1:length(timeMeteoTmp)
            % radiation
            tmp = squeeze(radiation(:,:,hourId))';
            meteoRadiationTmp(hourId,:) = interp2(x, y, tmp, siteLon, siteLat);
            
            % temperature
            tmp = squeeze(temperature(:,:,hourId))';
            meteoTemperatureTmp(hourId,:) = interp2(x, y, tmp, siteLon, siteLat);
            
            % precipitation
            tmp = squeeze(precipitation(:,:,hourId))';
            meteoPrecipitationTmp(hourId,:) = interp2(x, y, tmp, siteLon, siteLat);
            
            % relative humidity
            tmp = squeeze(humidity(:,:,hourId))';
            meteoHumidityTmp(hourId,:) = interp2(x, y, tmp, siteLon, siteLat);
            
            % speed
            tmp = squeeze(speed(:,:,hourId))';
            meteoSpeedTmp(hourId,:) = interp2(x, y, tmp, siteLon, siteLat);
        end
        
        meteoRadiation = [meteoRadiation; meteoRadiationTmp];
        meteoTemperature = [meteoTemperature; meteoTemperatureTmp];
        meteoPrecipitation = [meteoPrecipitation; meteoPrecipitationTmp];
        meteoHumidity = [meteoHumidity; meteoHumidityTmp];
        meteoSpeed = [meteoSpeed; meteoSpeedTmp];
    end
    meteoRadiation = interp1(timeMeteo, meteoRadiation, timeMeas);
    meteoTemperature = interp1(timeMeteo, meteoTemperature, timeMeas);
    meteoPrecipitation = interp1(timeMeteo, meteoPrecipitation, timeMeas);
    meteoHumidity = interp1(timeMeteo, meteoHumidity, timeMeas);
    meteoSpeed = interp1(timeMeteo, meteoSpeed, timeMeas);
else
    meteoFolder = 'meteo';
    meteoFile = [num2str(year) '.nc'];%SingleLevel_20190331
    
    filename = fullfile(dataFolder, meteoFolder, meteoFile);
    
    lon = ncread(filename, 'longitude');
    lat = ncread(filename, 'latitude');
    [x,y] = meshgrid(lon,lat);
    [xMeteo,yMeteo] = meshgrid(lon,lat);
    
    inputFormat = 'dd.MM.yyyy HH:mm';
    tmp = datetime('01.01.1900 00:00','InputFormat',inputFormat, 'TimeZone', 'UTC');
    formatOut = 'yyyymmddHH';
    currentTime = ncread(filename, 'time');
    meteoTime = tmp+double(currentTime)/24;
    
    timeMeteoTmp = hours(meteoTime - pollData.Date_time(1));
    timeMeteo = [timeMeteo; timeMeteoTmp];
    
    radiation =  ncread(filename, 'ssrd');
    radiation = max(diff(radiation,1,3), 0)/3600;
    radiation = cat(3, zeros(size(radiation,1),size(radiation,2)), radiation);
    radiation = squeeze(radiation);
    
    
    temperature = ncread(filename, 't2m')-273.15; % to celsius degree
    temperature =squeeze(temperature);
    
    
    precipitation = ncread(filename, 'tp')*1000;
    precipitation = max(diff(precipitation,1,3), 0);
    precipitation = cat(3, zeros(size(precipitation,1),size(precipitation,2)), precipitation);
    precipitation = squeeze(precipitation);
    
    dew = ncread(filename, 'd2m')-273.15; % to celsius degree
    humidity = 100*exp(17.625*dew./(243.04+dew))./exp((17.625*temperature)./(243.04+temperature));
    humidity = squeeze(humidity);
    
    u = ncread(filename, 'u10');
    v = ncread(filename, 'v10');
    speed = sqrt(u.^2+v.^2);
    speed = squeeze(speed);
end

%% road traffic
roadFile = fullfile(dataFolder, 'roadData', 'trafficVol.mat');
load(roadFile);
if(~predictionFlag)
    roadIntensity = zeros(size(newTable));
    for siteId = 1:size(roadIntensity,2)
        roadIntensity(:, siteId) = interp2(lonRoadGrid, latRoadGrid, trafficVol, sites{siteId, 4}, sites{siteId, 3});
    end
else
    roadIntensity = zeros(size(xnew));
    delteRoad = lonRoad(2)-lonRoad(1);
    for i=1:size(xnew,1)
        for j = 1:size(xnew,2)
            xtmp = xnew(i,j);
            ytmp = ynew(i,j);
            idl = round(max(((xtmp-Delta_x/2)-lonRoad(1)+delteRoad/2)/delteRoad+1,1));
            idr =round( min(((xtmp+Delta_x/2)-lonRoad(1)+delteRoad/2)/delteRoad+1,length(lonRoad)));
            
            idu = round(max((latRoad(1)+delteRoad/2-(ytmp+Delta_y/2))/delteRoad+1,1));
            idd = round(min((latRoad(1)+delteRoad/2-(ytmp-Delta_y/2))/delteRoad+1,length(latRoad)));
            
            %             idValid= find(lonRoadGrid<=(xtmp+Delta_x/2)& lonRoadGrid>=(xtmp-Delta_x/2)&latRoadGrid<=(ytmp+Delta_y)&latRoadGrid>=(ytmp-Delta_y));
            vtmp = trafficVol(idu:idd, idl:idr);
            roadIntensity(i,j) = sum(vtmp(:))/length(vtmp(:));
        end
    end
    %     roadIntensity = interp2(lonRoadGrid, latRoadGrid, trafficVol, xnew, ynew);
end


%% read population
popFilename = fullfile(dataFolder, 'pop', 'pop.tif');

[pop, popCoord] = geotiffread(popFilename);
cellLon = popCoord.CellExtentInLongitude;
cellLat = popCoord.CellExtentInLatitude;
lonPop = (popCoord.LongitudeLimits(1)+cellLon/2):cellLon:popCoord.LongitudeLimits(2);
latPop = (popCoord.LatitudeLimits(2)-cellLat/2):-cellLat:popCoord.LatitudeLimits(1);

if(~predictionFlag)
    popDensity = zeros(size(newTable));
    for siteId = 1:size(popDensity,2)
        popDensity(:, siteId) = interp2(lonPop, latPop, pop, sites{siteId, 4}, sites{siteId, 3});
    end
else
    popDensity = interp2(lonPop, latPop, pop, xnew, ynew);
end



if(~predictionFlag)
    %% combine all the data
    allData = table();
    % 1st cams
    allData.cams = camsNOX(:);
    
    % 2nd radiation
    allData.radiation = meteoRadiation(:);
    
    % 3rd temperature
    allData.temperature = meteoTemperature(:);
    
    % 4th precipitation
    allData.precipitation = meteoPrecipitation(:);
    
    % 5th humidty
    allData.humidity = meteoHumidity(:);
    
    % 6th speed
    allData.Speed = meteoSpeed(:);
    
    % 7th road
    allData.road = roadIntensity(:);
    
    % 8th pop
    % allData.pop = popDensity(:);
    
    % 9, 10, 11 time properties
    allData.hour = repmat(hour(pollData.Date_time),size(sites,1),1);
    allData.month = repmat(month(pollData.Date_time),size(sites,1),1);
    allData.weekday = repmat(weekday(pollData.Date_time),size(sites,1),1);
    
    % measurements of pollutant
    tmp = table2array(newTable);
    allData.measurements = tmp(:)./allData.cams(:);
    
    writetable(allData, fullfile(folder, [pollutant '_trainData.csv']));
else
    for hourId = 0:8759
        disp(hourId)
        %% combine all the data
        allData = table();
        % 1st cams
        id = find(timeCAMS == hourId);
        camsNOX = interp2(xCAMS, yCAMS, squeeze(camsAllData(:,:,id))', xnew, ynew);
        allData.cams = camsNOX(:);
        
        % 2nd radiation
        idMeteo = find(timeMeteo == hourId);
        meteoRadiation = interp2(xMeteo, yMeteo, squeeze(radiation(:,:,idMeteo))', xnew, ynew);
        allData.radiation = meteoRadiation(:);
        
        % 3rd temperature
        meteoTemperature = interp2(xMeteo, yMeteo, squeeze(temperature(:,:,idMeteo))', xnew, ynew);
        allData.temperature = meteoTemperature(:);
        
        % 4th precipitation
        meteoPrecipitation =  interp2(xMeteo, yMeteo, squeeze(precipitation(:,:,idMeteo))', xnew, ynew);
        allData.precipitation = meteoPrecipitation(:);
        
        % 5th humidty
        meteoHumidity = interp2(xMeteo, yMeteo, squeeze(humidity(:,:,idMeteo))', xnew, ynew);
        allData.humidity = meteoHumidity(:);
        
        % 6th speed
        meteoSpeed = interp2(xMeteo, yMeteo, squeeze(speed(:,:,idMeteo))', xnew, ynew);
        allData.Speed = meteoSpeed(:);
        
        % 7th road
        allData.road = roadIntensity(:);
        
        nSize = length(roadIntensity(:));
        % 8th pop
        % allData.pop = popDensity(:);
        ttmp = meteoTime(idMeteo);
        ttmp.TimeZone = 'Europe/Berlin';
        allData.hour = repmat(hour(ttmp),nSize,1);
        allData.month = repmat(month(ttmp),nSize,1);
        allData.weekday = repmat(weekday(ttmp),nSize,1);
        
        %     allData(isnan(allData)) = 0;
        %         writetable(allData, fullfile(folder, num2str(year), [num2str(hourId) '_predData.csv']));
        fileID = fopen(fullfile(folder, num2str(year), [num2str(hourId) '_' pollutantType '_predData.bin']),'w');
        fwrite(fileID,table2array(allData),'float');
        fclose(fileID);
    end
end
