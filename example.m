%% Example 1, plot time series of stuff
clear

time1 = 130; % ka 1950
time2 = -10; % ka 1950

% load orbital parameters
[tka, ecc, obl, lpe, pre] = getlaskar2004(5,'slice',[time1 time2]);
% [tka, ecc] = getlaskar2010(1,'slice',[time1 time2]); % Laskar2010 eccentricity solution only makes a difference for >20 Ma

% Calculate 65N summer
[inso] = insolationwm2(65, 90, 1367, ecc, obl, lpe);

% plot
figure(1)
clf

subplot(6,1,1)
plot(tka,rad2deg(obl))
ylabel(['Obliquity ( ',char(0176),' )'])
set(gca,'xdir','reverse')
set(gca,'xgrid','on')

subplot(6,1,2)
plot(tka,ecc)
ylabel(['Eccentricity'])
set(gca,'xdir','reverse')
set(gca,'xgrid','on')

subplot(6,1,3)
plot(tka,pre)
hold on
ylabel(['Precession index'])
set(gca,'xdir','reverse')
set(gca,'xgrid','on')
set(gca,'ygrid','on')

subplot(6,1,4)
plot(tka,inso)
set(gca,'xdir','reverse')
ylabel(['Mean daily summer solstice',newline,'insolation at 65N (W m^-^2)'])
set(gca,'xgrid','on')
set(gca,'ygrid','on')

% total J/m2 received at 65N for all days in the year
% where mean daily insolation is >= 275 W/m2
% (remake of Huybers et al (2006) figure 2C & E, to check if everything works)
subplot(6,1,5)
threshold = 275;
meaninso = NaN(size(tka));
intinso = NaN(size(tka));
days = 0:0.01:365.23; % because 0 = 365.24
for i = 1:numel(tka)
	sunlons = sday2sunlon(days,ecc(i),lpe(i),365.24); % solar longitude associated with each day
	insos = insolationwm2(65, sunlons, 1367, ecc(i), obl(i), lpe(i));
	meaninso(i) = mean(insos(insos>=threshold));
	secs = ( max(days(insos>=threshold)) - min(days(insos>=threshold)) ) * (24*60*60); % this approach works at high latitudes only
	intinso(i) = meaninso(i) * secs; % W/m2 -> J/m2
end
plot(tka,intinso/10^9)
xlabel('Age (ka b1950)')
ylabel('GJm^-^2 65N')
set(gca,'xdir','reverse')
grid on

% do same for 65S
subplot(6,1,6)
threshold = 275;
meaninso = NaN(size(tka));
intinso = NaN(size(tka));
days = 0:0.01:365.23; % because 0 = 365.24
for i = 1:numel(tka)
	sunlons = sday2sunlon(days,ecc(i),lpe(i),365.24);
	insos = insolationwm2(-65, sunlons, 1367, ecc(i), obl(i), lpe(i));
	meaninso(i) = mean(insos(insos>=threshold));
	secs = ( max(days(insos>=threshold)) - min(days(insos>=threshold)) ) * (24*60*60);
	intinso(i) = meaninso(i) * secs; % W/m2 -> J/m2
end
plot(tka,intinso/10^9)
xlabel('Age (ka b1950)')
ylabel('GJm^-^2 65S')
set(gca,'xdir','reverse')
grid on


%% Example 2: divergent heatmap comparing LGM Earth with 2000 AD Earth
clear

% grid
lats = [-89:1:89];
days = [0:0.1:365.2];

% 2000 AD Earth
[tka ecc obl lpe pre] = getlaskar2004(5,'slice',[-0.05 -0.05]); %#ok<*ASGLU> 2000 AD (-0.05 ka 1950)
sunlons = sday2sunlon(days,ecc,lpe,365.2);
insomov = NaN(numel(lats),numel(days));
lod1 = NaN(size(insomov));
for i = 1:numel(lats)
	[insomov(i,:) lod1(i,:)] = insolationwm2(lats(i), sunlons, 1367, ecc, obl, lpe);
end

% LGM Earth
[tka ecc obl lpe pre] = getlaskar2004(5,'slice',[21.950 21.950]); %#ok<*ASGLU> closest to 22 ka b1950 in Laskar data
sunlons = sday2sunlon(days,ecc,lpe,365.2);
insomesh2 = NaN(numel(lats),numel(days));
lod2 = NaN(size(insomov));
for i = 1:numel(lats)
	[insomesh2(i,:) lod2(i,:)] = insolationwm2(lats(i), sunlons, 1367, ecc, obl, lpe);
end

insomesh = insomov-insomesh2;
daymesh = lod1-lod2;

figure(2)
clf
plh = imagesc(days,lats,insomesh);
%[plc plh] = contourf(sunlon,lats,insomesh);
%hcl = clabel(plc,plh,'FontSize',8,'Color','w','Rotation',0);
hold on
grid on
plot([0 359],[0 0],'k:')
set(gca,'ydir','normal')
xticks([0 365.2/4 365.2/2 365.2*0.75 365.2]);
%xticklabels({'NH Spring','NH Summer','NH Autumn','NH Winter','NH Spring'})
cb = colorbar;
set(cb,'position',[0.93 0.11 0.02 .5])
title('Modern minus LGM W/m2')
%title('Hours of daylight by latitude and time of year')
%title(cb,'Hours')

xlabel('Day of year (relative to NH spring equinox)')
ylabel('Latitude (degrees N)')

%diverging colourmap (found on matlab central, lost link)
L = numel(insomesh);
indexValue = 0;     % value for which to set a particular color
topColor = [0.706, 0.016, 0.150];         % color for maximum data value (red = [1 0 0])
indexColor = [1 1 1];       % color for center colour in diverging colour map
bottomcolor = [0.230, 0.299, 0.754];      % color for minimum data value (blue = [0 0 1])
% Calculate where proportionally indexValue lies between minimum and
% maximum values
largest = max(max(insomesh));
smallest = min(min(insomesh));
index = L*abs(indexValue-smallest)/(largest-smallest);
% Create color map ranging from bottom color to index color
% Multipling number of points by 100 adds more resolution
customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
	linspace(bottomcolor(2),indexColor(2),100*index)',...
	linspace(bottomcolor(3),indexColor(3),100*index)'];
% Create color map ranging from index color to top color
% Multipling number of points by 100 adds more resolution
customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
	linspace(indexColor(2),topColor(2),100*(L-index))',...
	linspace(indexColor(3),topColor(3),100*(L-index))'];
customCMap = [customCMap1;customCMap2];  % Combine colormaps
colormap(customCMap)
clear customCMap1 customCMap2

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

clear customCMap

%% Example 3: megamovie (also as animated GIF)

clear

time1 = 30; % ka 1950
time2 = -10; % ka 1950

% load orbital parameters
[tka, ecc, obl, lpe, pre] = getlaskar2004(5,'slice',[time1 time2]);

% % make hires (comment out to not do it)
% tkahi = 1600:-1:0;
% ecc = interp1(tka,ecc,tkahi);
% obl = interp1(tka,obl,tkahi);
% lpe = interp1(tka,lpe,tkahi);
% pre = interp1(tka,pre,tkahi);
% tka = tkahi;

% grid
latres = 0.5;
lats = [89.75:-latres:-89.75]; % centre of pixel coordinates
dayres = 0.2;
days = [0.1:dayres:365.10]; % centre of pixel coordinates
areas = NaN(numel(lats),1);
for i = 1:numel(lats)
	areas(i,1) = areaquad(lats(i)-latres/2,0,lats(i)+latres/2,360,referenceEllipsoid('wgs84','meters'));
end

% make movie frames
insomov = NaN(numel(lats),numel(days),numel(tka));
daymov = NaN(size(insomov));
summerdn = sunlon2sday(90,ecc,lpe,365.2);
summerds = sunlon2sday(270,ecc,lpe,365.2);
%latsmesh = NaN(size(insomov));
for j = 1:numel(tka)
	disp([num2str(tka(j)), ' ka'])
	sunlons = sday2sunlon(days,ecc(j),lpe(j),365.2);
	for i = 1:numel(lats)
		[insomov(i,:,j) daymov(i,:,j)] = insolationwm2(lats(i), sunlons, 1367, ecc(j), obl(j), lpe(j));
		%latsmesh(i,:,j) = lats(i); % sanity check
	end
end

maxinso = max(insomov,[],'all');

domovie = 1; % enable to create the animated GIF
if domovie == 1
	figure(3)
	count = 0;
	for i = 1:size(insomov,3)
		clf
		plh = imagesc(days,lats,insomov(:,:,i));
		%[plc plh] = contourf(sunlon,lats,insomesh);
		%hcl = clabel(plc,plh,'FontSize',8,'Color','w','Rotation',0);
		hold on
		grid on
		plot([0 359],[0 0],'k:')
		set(gca,'ydir','normal')
		xticks([0 365.2/4 365.2/2 365.2*0.75 365.2]);
		%xticklabels({'NH Spring','NH Summer','NH Autumn','NH Winter','NH Spring'})
		cb = colorbar;
		set(cb,'position',[0.93 0.11 0.02 .5])
		title([num2str(round(tka(i))),' ka, Wm^-^2'])
		%title('Hours of daylight by latitude and time of year')
		%title(cb,'Hours')
		caxis([275 600]);
		
		plot([summerdn(i)],[65],'w+')
		plot([summerds(i)],[-65],'w+')
		
		xlabel('Day of year (relative to NH spring equinox)')
		ylabel('Latitude (degrees N)')
		
		% set figure size (cm)
		xSize = 20;
		ySize = 15;
		% set paper size (cm)
		set(gcf,'PaperUnits','centimeters')
		Y = xSize+1;
		X = ySize+1;
		set(gcf, 'PaperSize',[X Y])
		
		% put figure in centre of paper
		xLeft = (X-xSize)/2;
		yBottom = (Y-ySize)/2;
		set(gcf,'PaperPosition',[xLeft yBottom xSize ySize])
		% make background white
		set(gcf,'InvertHardcopy','on');
		set(gcf,'color',[1 1 1]);
		
		filename = 'deglaciation.gif';
		frame = getframe(gcf);
		im = frame2im(frame);
		[imind,cm] = rgb2ind(im,256);
		
		count = count + 1;
		if count == 1
			imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
		else
			imwrite(imind,cm,filename,'gif','WriteMode','append');
		end
		
	end
end

save('/media/stuff/orbits/insomov.mat','-v7.3','-nocompression','insomov','areas','dayres','tka','ecc','lpe','pre','lats');
% compression takes too long


