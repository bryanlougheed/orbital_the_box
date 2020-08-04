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
% con = Solar constant. Single numerical value or 1D array, w/m2. Leave empty, i.e. [], for 1367.
% ecc = Eccentricity. Numerical value(s). 1D array.
% obl = Obliquity. Numerical value(s), radians. 1D array.
% lpe = longitude of perihelion from moving equinox.
%		Numerical value(s), radians. 1D array.
%
% ecc, obl and lpe must be same dimensions.
%
% Output
% ======
%
% inso = Calculated insolation at top of atmosphere. W/m2
%        Array same size as ecc, obl and lpe.
%
% dayhrs = Hours of daylight.
%	   Array same size as ecc, obl and lpe.
%
% B.C. Lougheed, May 2020
% Matlab 2019a
%
% -----------------------------------------
% Insolation script unceremoniously stolen from 
% http://eisenman.ucsd.edu/code/daily_insolation.m, 
% modified to accept Laskar solutions and rewritten
% slightly in a way that I understand better.
%
% I added daylight hours following sunrise equation: 
% https://en.wikipedia.org/wiki/Sunrise_equation

if isempty(con) == 1
	con = 1367;
end

sunlon = deg2rad(sunlon);
lat = deg2rad(lat);

% declination angle of the sun
% https://en.wikipedia.org/wiki/Position_of_the_Sun
dsun = asin(sin(obl) .* sin(sunlon));

% Hour angle at sunrise/sunset
% https://en.wikipedia.org/wiki/Sunrise_equation
hangle = acos(-tan(lat).*tan(dsun));

% polar night / day. Berger (1978), eqquations a 8 and 9
hangle( abs(lat) >= pi/2 - abs(dsun) &  lat.*dsun  > 0  ) = pi;
hangle( abs(lat) >= pi/2 - abs(dsun) &  lat.*dsun <= 0  ) = 0;

% change lpe from heliocentric to geocentric: Appendix B in Berger et al. (2010), doi:10.1016/j.quascirev.2010.05.007
omega = lpe + pi; % add 180
omega(omega>2*pi) = omega(omega>2*pi) - 2*pi; % subtract 360 from stuff greater than 360

% Insolation: Berger 1978 eq (10)
inso  = con/pi.*(1+ecc.*cos(sunlon-omega)).^2 ./ ...
        (1-ecc.^2).^2 .* ...
        ( hangle.*sin(lat).*sin(dsun) + cos(lat).*cos(dsun).*sin(hangle) );
	
% Hours of daylight (https://en.wikipedia.org/wiki/Sunrise_equation)
dayhrs = abs(hangle - hangle*-1) / (2*pi/24);


