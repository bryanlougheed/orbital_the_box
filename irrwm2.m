function  [irr, dayhrs, tsi, rx] = irrwm2(lat, sunlon, con, ecc, obl, lpe, earthshape)
% [irr, dayhrs, tsi, rx] = irrwm2(lat, sunlon, con, ecc, obl, lpe, earthshape)
%
% Calculate daily average irradiance (W/m2) at top of atmosphere
% And also length of daytime (in hours)
%
% Input
% =====
%
% lat = Latitude (in degrees N, so minus for S) on Earth. 1D array.       _
% sunlon = Geocentric solar longitude (lambda). 1D array. (degrees, 0-359.9)
%          Northern spring equinox = 0
%          Northern summer solstice = 90
%          Northern autumn equinox = 180
%          Northern winter solstice = 270
% con = Solar constant. Single numerical value or 1D array, W/m2. Leave empty, i.e. [], for 1367.
% ecc = Eccentricity. Numerical value(s). 1D array.
% obl = Obliquity. Numerical value(s), radians. 1D array.
% lpe = longitude of perihelion from moving equinox. (omega-bar)
%       Numerical value(s), radians. 1D array.
% earthshape = Optional. Shape of Earth, enter string 'sphere' or 'wgs84' (optional, default is 'sphere')
%
% ecc, obl and lpe must be same dimensions.
%
% Output
% ======
%
% irr    = Calculated mean daily (24 hr) irradiance at top of atmosphere. W/m2
%          Array same size as ecc, obl and lpe.
%
% dayhrs = Hours of daylight.
%          Array same size as ecc, obl and lpe.
%
% tsi  =   Calculated mean daily irrlation at top of atmosphere assuming 90 degree angle of incidence. W/m2
%          Insensitive to latitude, obliquity or earthshape.
%          Array same size as ecc, obl and lpe. 
%
% rx     = Distance from Sun (AU). Insensitive to latitude, obliquity or earthshape.
%
% -----------------------------------------
%
% B.C. Lougheed, May 2020, Matlab 2019a
% Updated with wgs84 oblateness Sep. 2020
%
% -----------------------------------------
% irr in Wm2 based on equations in Berger (1978).
%
% Part of script (specifically polar night and day) uses Ian Eisenman's
% Matlabification of Berger (1978) equations 8 and 9 (see comments in script)
% taken from here: http://eisenman.ucsd.edu/code/daily_insolation.m
%
% I added ability to take oblateness of Earth into account,
% validated using Van Hemelrijck (1983) solution. (see comments in script)
%
% I added daylight hours output following sunrise equation: 
% https://en.wikipedia.org/wiki/Sunrise_equation
%
% I added irr90 by calculating distance from sun following Meeus (1998).
% Meeus, J., (1998). Astronomical Algorithms, 2nd ed. Willmann-Bell, Inc., Richmond, Virginia. (specifically Chapter 30).


if isempty(con) == 1
	con = 1367;
end

if nargin < 7
	earthshape = 'sphere';
end

sunlon = deg2rad(sunlon);

if strcmpi(earthshape,'sphere') == 1
	lat = deg2rad(lat);
elseif strcmpi(earthshape,'wgs84') == 1 % other planets could be done similarly, I assume
	f = 1/298.257223563; % wgs84 flattening value
	re = 6378137; % wgs84 equatorial radius (metres)
	rp = re*(1-f); % calculate polar radius
	% replace geocentric lat with geographic lat in all subsequent calculations
	lat = atan( (re/rp)^2 * tan(deg2rad(lat))); 
	% note: It is also possible to use the very long equation 11 in Van Hemelrijck (1983) to calculate
	% irradiance for an oblate Earth using an angle correction upon the irradiance calculated using
	% geocentric latitude. I did that also and found that it gave exactly the same outcome as simply
	% subsituting geographic latitude directly into the Berger (1978) equation 10 below. Since the latter way
	% is computationally more efficient, I use that.
end
	
% declination angle of the sun
% https://en.wikipedia.org/wiki/Position_of_the_Sun
dsun = asin(sin(obl) .* sin(sunlon));

% Hour angle at sunrise/sunset
% https://en.wikipedia.org/wiki/Sunrise_equation
hangle = acos(-tan(lat).*tan(dsun));
% polar night / day. Berger (1978), eq. 8 and 9
% following two lines done following the Ian Eisenman script:
hangle( abs(lat) >= pi/2 - abs(dsun) &  lat.*dsun  > 0  ) = pi; % polar day
hangle( abs(lat) >= pi/2 - abs(dsun) &  lat.*dsun <= 0  ) = 0;  % polar night

% Hours of daylight (https://en.wikipedia.org/wiki/Sunrise_equation)
dayhrs = abs(hangle - hangle*-1) / (2*pi/24); % assuming exactly 24 hours in a day

% change lpe from heliocentric to geocentric (i.e. change from omega-bar to omega): 
omega = lpe + pi; % add 180 degrees, see e.g. Appendix B in Berger et al. (2010), doi:10.1016/j.quascirev.2010.05.007
omega(omega>=2*pi) = omega(omega>=2*pi) - 2*pi; % subtract 360 from stuff >= 360 (put back in same orbital period)

% Irradiance: Berger (1978) eq (10)
irr  = con/pi.*(1+ecc.*cos(sunlon-omega)).^2 ./ (1-ecc.^2).^2 .* ( hangle.*sin(lat).*sin(dsun) + cos(lat).*cos(dsun).*sin(hangle) );

% % Van Hemelrijck (1983) fancier method for oblate Earth
% % in this case, lat is the geocentric latitude, while using hangle that was calculated using geographic latitude
% vangle = atan((1-f)^-2 * tan(lat))-lat; % Van Hemelrijck (1983) eq. 9, f is wgs84 flattening
% % Van Hemelrijck (1983) eq. 11 second term	
% hemelterm = (cos(vangle).*(hangle.*sin(lat).*sin(dsun)+sin(hangle).*cos(lat).*cos(dsun))+sin(vangle).*(-tan(lat).*(hangle.*sin(lat).*sin(dsun)+sin(hangle).*cos(lat).*cos(dsun))+hangle.*sin(dsun).*sec(lat)));
% % irradiation: Berger (1978) eq. 10, but replace final term with Van Hemelrijck (1983) eq. 11 second term.
% irr  = con/pi.*(1+ecc.*cos(sunlon-omega)).^2 ./ (1-ecc.^2).^2 .* hemelterm;

% Calculate rx and tsi
% First, calculate distance from Sun in AU
veq = 2*pi - omega; % v (true anomaly) of spring equinox relative to perihelion
vx = veq + sunlon; % v (true anomaly) of inputted sunlon relative to perihelion
vx(vx>2*pi) = vx(vx>2*pi) - 2*pi; % put back in 0-360 range 
%rp = (1-ecc); % perihelion distance in AU
%ra = (1+ecc); % aphelion distance in AU
%a = (rp+ra)/2; % semi-major axis in AU (always 1, duh)
rx = 1*(1-ecc.^2) ./ (1+ecc.*cos(vx)); % Eq. 30.3 in Meeus (1998)
% calculate tsi as function of con relative to 1 AU
tsi = con * (1./rx).^2;

end