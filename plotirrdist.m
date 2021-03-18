function plotirrdist(obl,ecc,lpe,con,earthshape,savename)
% plotirrdist(obl,ecc,lpe,con,earthshape,savename)
%
% Makes a plot of distribution of Earth's irradiance throughout the year
% and saves copy to hard drive.
%
% Input
% -----
% obl = obliquity in radians
% ecc = eccentricity of the ellipse
% lpe = longitude of perihelion as given by, e.g., Laskar: 
%       omega-bar (i.e. relative to NH autumn equinox) in radians.
% con = solar constant in W/m2 (e.g. 1361)
% earthshape = 'wgs84' or 'sphere'
% savename = filename of the png file to saved (string)
%
%
% Output
% ------
% A plot of irradiance profile of the Earth in the current figure window,
% or new figure window, if none is open.
%
% B.C. Lougheed, December 2020

latres = 0.5;
dayres = 0.2;
ndays = 365.2;

lats = [90-(latres/2):-latres:-90+(latres/2)];
days = 0:dayres:ndays-dayres;
irrmat = NaN(numel(lats),numel(days));
daymat = NaN(numel(lats),numel(days));

% prep irradiance matrix
sunlons = sday2sunlon(days,ecc,lpe,ndays);
soleqdays = sunlon2sday([90 180 270],ecc,lpe,ndays);
for i = 1:numel(lats)
	[irrmat(i,:), daymat(i,:), ~, rx] = irrwm2(lats(i), sunlons, con, ecc, obl, lpe, earthshape);
end
annmeans = mean(irrmat,2);


% plot the irradiance
clf
axirr = axes('position',[0.05 0.14 0.75 0.54]);
contourf(days,lats,irrmat,[0:25:650],'edgecolor','none');
hold on
caxis([0 650])
hcb = colorbar;
set(hcb,'location','southoutside')
set(hcb,'position',[0.05 0.055 0.4 0.02])
htb = annotation('textbox',[0.04 0.014 0.4 0.02],'string',['Daily (24 hr) Q_{mean} (W m^-^2)']);
set(htb,'linestyle','none')
colormap(parula(numel(0:25:650)-1))
% plot the days of the solstices and equinoxes
for i = 1:3
	plot([soleqdays(i) soleqdays(i)],ylim,'k--','color',[0.3 0.3 0.3])
end
% set axis stuff
set(gca,'tickdir','out')
set(gca,'ydir','normal')
ylabel(['Latitude (',char(0176),'N)'])
xlabel('Day of year since northern hemisphere spring equinox')
ylim([-90 90])
yticks([-90:30:90])
box off
xlims = xlim;

% plot the distance from Sun
axdis = axes('position',[0.05 0.76 0.75 0.1970]);
hpdis = plot(days,rx*1.495978707*10^11/1000/10^6);
hold on
xlim(xlims)
ylim([140 160])
for i = 1:3
	plot([soleqdays(i) soleqdays(i)],ylim,'k--','color',[0.3 0.3 0.3])
end
%xticklabels([])
yticks([140 150 160])
ylabel('Distance from Sun (10^6 km)')
xlabel('Day of year since northern hemisphere spring equinox')
set(gca,'xaxislocation','top')
set(gca,'yaxislocation','left')
box off
set(axdis,'ycolor',hpdis.Color);

% plot orbital speed
axvel = axes('position',[0.05 0.76 0.75 0.1970]);
plot(nan,nan)
hold on
hpvel = plot(days, earthspeed(rx)/1000 );
xlim(xlims)
ylim([28 32])
%for j = 1:3
%	plot([soldays(i,j) soldays(i,j)],ylim,'k--')
%end
xticklabels([])
yticks([28 30 32])
ylabel(['Orbital speed',newline,'(km s^-^1)'])
set(gca,'yaxislocation','right')
set(gca,'color','none')
box off
set(axvel,'ycolor',hpvel.Color);

% plot average annual irradiance
axave = axes('position',[0.827 0.14 0.13 0.54]);
plot(annmeans,lats,'k-')
ylim([-90 90])
xlim([160 420])
xticks(linspace(160,420,3))
yticks([-90:30:90])
set(gca,'yaxislocation','right')
xlabel('Annual Q_{mean} (W m^-^2)')
ylabel(['Latitude (',char(0176),'N)'])

% season lengths
axsea = axes('position',[0.05 0.68 0.75 0.08]);
text(0,0.5,['Equinox',newline,'\lambda = 0',char(0176)],'horizontalalignment','center','verticalalignment','middle','color',[0.3 0.3 0.3])
text(soleqdays(1),0.5,['Solstice',newline,'\lambda = 90',char(0176)],'horizontalalignment','center','verticalalignment','middle','color',[0.3 0.3 0.3])
text(soleqdays(2),0.5,['Equinox',newline,'\lambda = 180',char(0176)],'horizontalalignment','center','verticalalignment','middle','color',[0.3 0.3 0.3])
text(soleqdays(3),0.5,['Solstice',newline,'\lambda = 270',char(0176)],'horizontalalignment','center','verticalalignment','middle','color',[0.3 0.3 0.3])
text(ndays,0.5,['Equinox',newline,'\lambda = 360',char(0176)],'horizontalalignment','center','verticalalignment','middle','color',[0.3 0.3 0.3])
text(mean([0 soleqdays(1)]),0.5,[num2str(soleqdays(1)-0,'%0.2f'),' days'],'horizontalalignment','center','color',[0.3 0.3 0.3])
text(mean([soleqdays(1) soleqdays(2)]),0.5,[num2str(soleqdays(2)-soleqdays(1),'%0.2f'),' days'],'horizontalalignment','center','color',[0.3 0.3 0.3])
text(mean([soleqdays(2) soleqdays(3)]),0.5,[num2str(soleqdays(3)-soleqdays(2),'%0.2f'),' days'],'horizontalalignment','center','color',[0.3 0.3 0.3])
text(mean([soleqdays(3) ndays]),0.5,[num2str(ndays-soleqdays(3),'%0.2f'),' days'],'horizontalalignment','center','color',[0.3 0.3 0.3])
xlim(xlims)
ylim([0 1])
axis off

% uniform fonts
set(findall(gcf,'-property','FontSize'),'FontSize',7)

% save to file
print2png(gcf, savename, 19, 15, 150)

	% nested function
	function print2png(fig, filename, X, Y, dpi, xmarg, ymarg)
		% print2png(fig, filename, X, Y, dpi, xmarg, ymarg)
		%
		% Required:
		% ---------
		% fig = figure number of figure object, e.g. gcf
		% filename = string
		% X = PNG width (cm)
		% Y = PNG height (cm)
		% dpi = resolution in dots per inch
		%
		% Optional (will be zero if not given):
		% -------------------------------------
		% xmarg = sum of horizontal margins on both sides (cm)
		% ymarg = sum of vertical margins on both sides (cm)
		%
		% -----------------------
		% B.C. Lougheed, Oct 2020
		
		if nargin < 6
			xmarg = 0;
			ymarg = 0;
		end
		
		% set figure size (cm)
		xSize = X-xmarg;
		ySize = Y-ymarg;
		
		% set paper size (cm)
		set(gcf,'PaperUnits','centimeters')
		set(gcf, 'PaperSize',[X Y])
		
		% put figure in centre of paper
		xLeft = (X-xSize)/2;
		yBottom = (Y-ySize)/2;
		set(gcf,'PaperPosition',[xLeft yBottom xSize ySize])
		
		% make background white
		set(gcf,'InvertHardcopy','on');
		set(gcf,'color',[1 1 1]);
		
		% save to hard drive
		print(figure(fig), '-dpng', ['-r',num2str(dpi)], [filename]);
		
	end


end


