function plotorbit(ecc,lpe,savename)
% plotorbit(ecc,lpe)
%
% Creates a scale plot of the Earth's orbit and saves to hard drive.
%
% Input
% -----
% ecc = eccentricity of the ellipse
% lpe = longitude of perihelion as given by, e.g., Laskar: 
%       omega-bar (i.e. relative to NH autumn equinox) in radians.
% savename = the name of the output .png file (string)
%
% Output
% ------
% A scale plot of the orbit in current plotting window, or new plotting
% window if there is none open.
%
% B.C. Lougheed, Jan. 2021

% true anomaly v of seasons
omega = lpe + pi; % perihelion
nse = 2*pi - omega; % nh spring equinox
nss = nse + pi/2; % nh summer solstice
nae = nss + pi/2; % nh autumn equinox
nws = nae + pi/2; % nh winter solstice
% convert true anomaly v to ellipse angle E (Meeus 1998 Equation 30.1 solve for E)
nse = 2 * atan( tan(nse/2) .* sqrt((1-ecc)./(1+ecc)) );
nss = 2 * atan( tan(nss/2) .* sqrt((1-ecc)./(1+ecc)) );
nae = 2 * atan( tan(nae/2) .* sqrt((1-ecc)./(1+ecc)) );
nws = 2 * atan( tan(nws/2) .* sqrt((1-ecc)./(1+ecc)) );

% semi-minor and semi-major axis of ellipse
a = 149.5978707; % au in 10^6 km
b = a*(1-ecc^2)^0.5;
% perihelion and aphelion distances
rper = a * ( 1 - ecc );
raph = a * ( 1 + ecc );
fd = a-rper; % focal distance offset from centre of ellipse

% ellipse circumfrence points x and y coords
E = linspace(0, 2*pi, 1000); % central angle E of ellipse, 0 to 360
xell = a*cos(E)-fd; % -fd to place sun at centre of image
yell = b*sin(E);

% plot
clf

% plot the ellipse
plot(xell, yell,'k-','color',[0.6 0.6 0.6], 'LineWidth', 1);
hold on

% plot the sun
plot(0,0,'yo','markersize',15,'markerfacecolor',[255 221 66]/255,'markeredgecolor',[255 143 66]/255)

% plot the seasons
ofs = 1.13; % text  offset
% text rotation
nserot = -90; if nse < 0; nserot =  90; end
nssrot = -90; if nss < 0; nssrot =  90; end
naerot = -90; if nae < 0; naerot =  90; end
nwsrot = -90; if nws < 0; nwsrot =  90; end
% plot quarter orbit lines
plot([a*cos(nws),a*cos(nss)]-fd,[b*sin(nws),b*sin(nss)],'r:')
plot([a*cos(nse),a*cos(nae)]-fd,[b*sin(nse),b*sin(nae)],'r:')
% spring equinox
plot(a*cos(nse)-fd, b*sin(nse),'ko','markersize',7,'markerfacecolor',[0 82 162]/255,'markeredgecolor',[00 100 00]/255)
text(a*cos(nse)*ofs-fd, b*sin(nse)*ofs, ['NH spring equinox',newline','\lambda = 0',char(0176)],'horizontalalignment','center','rotation',rad2deg(nse)+nserot)
% summer solstice
plot(a*cos(nss)-fd, b*sin(nss),'ko','markersize',7,'markerfacecolor',[0 82 162]/255,'markeredgecolor',[00 100 00]/255)
text(a*cos(nss)*ofs-fd, b*sin(nss)*ofs, ['NH summer solstice',newline','\lambda = 90',char(0176)],'horizontalalignment','center','rotation',rad2deg(nss)+nssrot)
% autumn equinox
plot(a*cos(nae)-fd, b*sin(nae),'ko','markersize',7,'markerfacecolor',[0 82 162]/255,'markeredgecolor',[00 100 00]/255)
text(a*cos(nae)*ofs-fd, b*sin(nae)*ofs, ['NH autumn equinox',newline','\lambda = 180',char(0176)],'horizontalalignment','center','rotation',rad2deg(nae)+naerot)
% winter solstice
plot(a*cos(nws)-fd, b*sin(nws),'ko','markersize',7,'markerfacecolor',[0 82 162]/255,'markeredgecolor',[00 100 00]/255)
text(a*cos(nws)*ofs-fd, b*sin(nws)*ofs, ['NH winter solstice',newline','\lambda = 270',char(0176)],'horizontalalignment','center','rotation',rad2deg(nws)+nwsrot)

% plot the per and aph distance
% per
ar(1) = arrow([xell(1) 0],[xell(1)-rper 0]);
ar(2) = arrow([xell(1)-rper 0],[xell(1) 0]);
text(mean([xell(1),0]),13,['Perihelion',newline,num2str(rper,'%0.1f'),' million km'],'horizontalalignment','center')
% aph
ar(3) = arrow([a*cos(pi)-fd 0],[0 0]);
ar(4) = arrow([0 0],[a*cos(pi)-fd 0]);
text(mean([a*cos(pi)-fd,0]),13,['Aphelion',newline,num2str(raph,'%0.1f'),' million km'],'horizontalalignment','center')

xlim([-170 170])
ylim([-170 170])

axis off
set(findall(gcf,'-property','FontSize'),'FontSize',6)
print2png(gcf, savename, 7.5, 7.5, 150)

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