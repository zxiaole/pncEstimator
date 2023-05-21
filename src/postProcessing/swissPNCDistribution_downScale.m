%% Aug. 22, 2021
clear
close all

%%
year = 2019;
Nx = 54;
Ny = 22;
fontsize = 16;
diffFlag = 0;
charLength = 1/30;
pollutant = 'NO2';

%% load average
load(['pollutants/downScale/' pollutant '_' num2str(year) '_annualAverage.mat']);

[xnew, ynew] = meshgrid(lonNew, latNew);
Delta_x = mean(diff(lonNew));
Delta_y = abs(mean(diff(latNew)));
%%

figure('position', [100 100 800 600])
mapFile ='SwissGridData/gadm36_CHE_shp/gadm36_CHE_0.shp';

roi = shaperead(mapFile);

names = {roi.NAME_0};

% regionId = find(strcmp(names, 'Switzerland'));


for regionId=1:length(roi)
    rx = roi(regionId).X;
    ry = roi(regionId).Y;
    
    plot(rx, ry, 'k-', 'linewidth',1.5)
    
    blockNum = sum(isnan(rx));
    blockId = find(isnan(ry));
    blockId = [0 blockId];
    
    x_min = min(lonNew);
    y_min = min(latNew);
    Ny = size(avgConc,1);
    Nx = size(avgConc,2);
    maskTmp = 0;
    for i =1:blockNum
        % convert to image coordinates
        ix = double((rx((blockId(i)+1):(blockId(i+1)-1)) - x_min)/Delta_x + 1);
        iy = Ny- double((ry((blockId(i)+1):(blockId(i+1)-1)) - y_min)/Delta_y + 1);
        maskTmp = min(maskTmp + poly2mask(ix,iy,Ny,Nx),1);
        
    end
    
end


imagesc(lonNew, latNew, log10(avgConc).*maskTmp, 'Alphadata', maskTmp,[0 2])

daspect([1 cos(mean(latNew)/180*pi) 1])
axis xy

colormap(jet)
cb =  cbarrow2();

h=gcf;
c=get(h,'children'); % Find allchildren
cb=findobj(h,'Tag','Colorbar'); % Find thecolorbar children
cb.FontSize = fontsize;
barTicks = cb.Ticks;
for i=1:length(barTicks)
    cb.TickLabels{i} = ['10^{' num2str(barTicks(i)) '}'];
    cb.FontName = 'Arial';
end


% cbLabel = typeUnits{dataId};
cb.Label.String = {[pollutant ' concentration (μg\cdotm^{-3})']};

ylabel('Latitude (degree)');
xlabel('Longitude (degree)')

set(gca, 'FontName', 'Arial', ...
    'FontSize', fontsize);

set(gcf,'PaperPositionMode','auto')

print(['figures/DownScale_' pollutant],'-dpng','-r300')
