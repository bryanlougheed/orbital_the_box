function [wm2, days, dayhr] = wm2intraday(lat,ecc,obl,lpe,con,daysinyear,dayres)
% [wm2, days, dayhr] = wm2intraday(lat,ecc,obl,lpe,con,daysinyear,dayres)
%
% Calculate a tropical year's worth of intraday W/m2 for a particular 
% latitude and orbital configuration. Assumes that northern spring equinox 
% will occur at day 0.0 (i.e. exactly local midnight on the first day of the
% tropical year). Script takes equation of time into account.
%
% Input
% -----
% lat = latitude in degrees
% ecc = eccentricity of the ellipse
% obl = obliquity in radians
% lpe = longitude of perihelion as given by, e.g., Laskar et al.: 
%       omega-bar (i.e. relative to NH autumn equinox) in radians.
% con = solar constant in W/m2, leave empty for 1361 W/m2
% daysinyear = number of mean solar days in the year, leave emtpy for 365.240
% dayres = mean solar day resolution for analysis, leave empty for 0.001
%
% Output
% ------
% wm2 = vector of W/m2 for every day interval calculated
% days = all the solar day intervals. 0 corresponds to NH spring equinox.
% dayhr = hour of the mean solar day, same dimensions as wm2 and days
%
% B.C. Lougheed, April 2023, Matlab 2019a

if nargin < 5
	con = 1361;
	daysinyear = 365.240;
	dayres = 0.001;
end

if isempty(con); con = 1361; end
if isempty(daysinyear); daysinyear = 365.240; end
if isempty(dayres); dayres = 0.001; end

% create time series of mean solar day fractions
days = [0:dayres:daysinyear-dayres];
dayhr = (days - floor(days))*24;
		
% calculate Earth's solar longitude (i.e. lambda) and equation of time for each day fraction,
[sunlon eot] = time2sunlon(days,ecc,lpe,daysinyear,obl); % input radians, output degrees... it made sense at the time, sorry :)
sunlon = deg2rad(sunlon);

% get local apparent solar hr (correct for eot)
solhr = (eot/60) + dayhr;
solhr(solhr<0) = solhr(solhr<0) + 24;
solhr(solhr>24) = solhr(solhr>24) - 24;

% declination of the sun, https://en.wikipedia.org/wiki/Position_of_the_Sun
dsun = asin(sin(obl) .* sin(sunlon));

% local hour angle (-180 to +180 deg, midday = 0 deg)
hangles = (2*pi/24) * (solhr - 12);

% solar elevation, https://en.wikipedia.org/wiki/Solar_zenith_angle
elev = asin( sin(dsun).*sin(deg2rad(lat)) + cos(dsun).*cos(deg2rad(lat)).*cos(hangles) );

% Calculate rx (distance from sun) and corresponding tsi
% First, calculate distance from Sun in AU
omega = lpe + pi;
veq = 2*pi - omega; % v (true anomaly) of NH spring equinox relative to perihelion
vx = veq + sunlon; % v (true anomaly) of inputted sunlon relative to perihelion
vx(vx>2*pi) = vx(vx>2*pi) - 2*pi; % put back in 0-360 range
rx = 1*(1-ecc.^2) ./ (1+ecc.*cos(vx)); % Eq. 30.3 in Meeus (1998)
tsi = con * (1./rx).^2; % calculate tsi as function of con relative to 1 AU

% calc W/m2, vertical component of tsi
wm2 = tsi .* sin(elev);
wm2(wm2<0) = 0; % sun under horizon, night time
	
end


