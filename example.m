time1 = 800; % ka 1950
time2 = -200; % ka 1950

% load orbital parameters
[tka ecc obl lpe pre] = getlaskar2004(5,'slice',[time1 time2]);
%[tka ecc] = getlaskar2010(1,'slice',[time1 time2]); % Laskar2010 eccentricity solution only makes a difference for >20 Ma


% Calculate 65N summer
[inso1] = insolationwm2(65, 90, 1367, ecc, obl, lpe);
% Calculate 65N spring
[inso2] = insolationwm2(0, 90, 1367, ecc, obl, lpe);
% Calculate 65N winter
[inso3] = insolationwm2(65, 270, 1367, ecc, obl, lpe);
% Calculate 65S summer
[inso4] = insolationwm2(-65, 270, 1367, ecc, obl, lpe);

% plot
figure(54)
clf

subplot(4,1,1)
plot(tka,rad2deg(obl))
ylabel(['Obliquity ( ',char(0176),' )'])
set(gca,'xdir','reverse')
set(gca,'xgrid','on')

subplot(4,1,2)
plot(tka,ecc)
ylabel(['Eccentricity'])
set(gca,'xdir','reverse')
set(gca,'xgrid','on')


subplot(4,1,3)
plot(tka,pre)
hold on
%plot([min(tka) max(tka)],[0 0],'k-')
ylabel(['Precession index'])
set(gca,'xdir','reverse')
set(gca,'xgrid','on')
set(gca,'ygrid','on')


subplot(4,1,4)
plot(tka,inso1)
%plot(tka,inso2)
%plot(tka,inso3)
set(gca,'xdir','reverse')
ylabel(['Mean daily summer solstice',newline,'insolation at 65N (W m^-^2)'])
xlabel('Age (ka, thousands of years)')
set(gca,'xgrid','on')
set(gca,'ygrid','on')




%% make a heatmap of earth insolation for a particular year
clear

% grid
lats = [-89:1:89];
sunlon = [0:1:359];


% load orbital parameters
[tka ecc obl lpe pre] = getlaskar2004(5,'slice',[-0.05 -0.05]); %#ok<*ASGLU> 2000 AD (-0.05 ka 1950) is closest Laskar year
insomesh1 = NaN(numel(lats),numel(sunlon));
lod = NaN(size(insomesh1));
for i = 1:numel(lats)
	for j = 1:numel(sunlon)
		[insomesh1(i,j) lod(i,j)] = insolationwm2(lats(i), sunlon(j), 1367, ecc, obl, lpe);
	end
end

% load orbital parameters
[tka ecc obl lpe pre] = getlaskar2004(5,'slice',[21.950 21.950]);
insomesh2 = NaN(numel(lats),numel(sunlon));
for i = 1:numel(lats)
	for j = 1:numel(sunlon)
		[insomesh2(i,j)] = insolationwm2(lats(i), sunlon(j), 1367, ecc, obl, lpe);
	end
end

insomesh = insomesh1-insomesh2;

clf
plh = imagesc(sunlon,lats,insomesh);
%[plc plh] = contourf(sunlon,lats,insomesh);
%hcl = clabel(plc,plh,'FontSize',8,'Color','w','Rotation',0);
hold on
grid on
plot([0 359],[0 0],'k:')
set(gca,'ydir','normal')
xticks([0 90 180 270 359]);
xticklabels({'NH Spring','NH Summer','NH Autumn','NH Winter','NH Spring'})
cb = colorbar;
set(cb,'position',[0.93 0.11 0.02 .5])
title('Modern minus LGM W/m2')
%title('Hours of daylight by latitude and time of year')
%title(cb,'Hours')

xlabel('Time of year')
ylabel('Latitude (degrees)')

%diverging colormap
% data = insomesh; % create example data set with values ranging from 0 to 3.6
% L = numel(data);
% indexValue = 1;     % value for which to set a particular color
% topColor = [0.706, 0.016, 0.150];         % color for maximum data value (red = [1 0 0])
% indexColor = [0.99 0.99 0.99];       % color for indexed data value (white = [1 1 1])
% bottomcolor = [0.230, 0.299, 0.754];      % color for minimum data value (blue = [0 0 1])
% % Calculate where proportionally indexValue lies between minimum and
% % maximum values
% largest = max(max(data));
% smallest = min(min(data));
% index = L*abs(indexValue-smallest)/(largest-smallest);
% % Create color map ranging from bottom color to index color
% % Multipling number of points by 100 adds more resolution
% customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
%             linspace(bottomcolor(2),indexColor(2),100*index)',...
%             linspace(bottomcolor(3),indexColor(3),100*index)'];
% % Create color map ranging from index color to top color
% % Multipling number of points by 100 adds more resolution
% customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
%             linspace(indexColor(2),topColor(2),100*(L-index))',...
%             linspace(indexColor(3),topColor(3),100*(L-index))'];
% customCMap = [customCMap1;customCMap2];  % Combine colormaps
colormap(parula)
%colormap(flipud(brewer(24)))

%plot2raster(gca, {plh}, 'bottom', 300);

% set figure size (cm)
xSize = 20;
ySize = 10;
% set paper size (cm)
set(gcf,'PaperUnits','centimeters')
Y = ySize+0;
X = xSize+0;
set(gcf, 'PaperSize',[X Y])

% put figure in centre of paper
xLeft = (X-xSize)/2;
yBottom = (Y-ySize)/2;
set(gcf,'PaperPosition',[xLeft yBottom xSize ySize])
% make background white
set(gcf,'InvertHardcopy','on');
set(gcf,'color',[1 1 1]);

print(gcf, '-dpdf', '-painters', ['insotimemap.pdf']);
