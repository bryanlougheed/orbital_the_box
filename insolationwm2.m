function  [inso, dayhrs] = insolationwm2(lat, sunlon, con, ecc, obl, lpe)
% [inso dayhrs] = insolationwm2(lat, sunlon, con, ecc, obl, lpe)
%
% Calculate daily average insolation (W/m2) at top of atmosphere
% And also length of daytime (in hours)
%
% Input
% =====
%
% lat = Latitude (in degrees N) on Earth. 1D array.
% sunlon = Geocentric solar longitude of the sun. 1D array. (degrees, 0-359)
%		Northern spring equinox = 0
%		Northern summer solstice = 90
%		Northern autumn equinox = 180
%		Northern winter solstice = 270
% con = Solar constant. Single numerical value, w/m2. Leave empty, i.e. [], for 1367.
% ecc = Eccentricity. Numerical value(s). 1D array.
% obl = Obliquity. Numerical value(s), radians. 1D array.
% lpe = longitude of perihelion from moving equinox.
%		Numerical value(s), radians. 1D array.
%
% ecc, obl and per must be same dimensions.
%
% Output
% ======
%
% inso = Calculated insolation at top of atmosphere. W/m2
%        Array same size as ecc, obl and per.
%
% dayhrs = Hours of daylight.
%	   Array same size as ecc, obl and per.
%
% B.C. Lougheed, May 2020
% Matlab 2019a
%
% -----------------------------------------
% Insolation script unceremoniously stolen from 
% http://eisenman.ucsd.edu/code/daily_insolation.m
% and modified to accept Laskar solutions.
%
% Daylight hours added following Wikipedia: 
% https://en.wikipedia.org/wiki/Sunrise_equation

if isempty(con) == 1
	con = 1367;
end

sunlon = deg2rad(sunlon);
lat = deg2rad(lat);

% declination angle of the sun
dsun = asin(sin(obl) .* sin(sunlon));

% Hour angle at sunrise/sunset
hangle = acos(-tan(lat).*tan(dsun));

% polar night / day. Berger 1978 eqn (8),(9)
hangle( abs(lat) >= pi/2 - abs(dsun) &  lat.*dsun  > 0  ) = pi;
hangle( abs(lat) >= pi/2 - abs(dsun) &  lat.*dsun <= 0  ) = 0;

% change lpe from heliocentric to geocentric: Appendix B in Berger et al. (2010), doi:10.1016/j.quascirev.2010.05.007
glpe = lpe + pi; % add 180
glpe(glpe>2*pi) = glpe(glpe>2*pi) - 2*pi; % subtract 360 from stuff greater than 360

% Insolation: Berger 1978 eq (10)
inso  = con/pi*(1+ecc.*cos(sunlon-glpe)).^2 ./ ...
        (1-ecc.^2).^2 .* ...
        ( hangle.*sin(lat).*sin(dsun) + cos(lat).*cos(dsun).*sin(hangle) );
	
% Hours of daylight (https://en.wikipedia.org/wiki/Sunrise_equation)
dayhrs = abs(hangle - hangle*-1) / (2*pi/24);


