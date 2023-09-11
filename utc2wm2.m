function [wm2 sunlon eot lappst elev] = utc2wm2(dtime, lat, lon, con)
% [wm2 sunlon eot lappst elev] = utc2wm2(dtime, lat, lon, con)
%
% Get the top of atmosphere irradiance for a particular lat/lon 
% on Earth at a precise UTC datetime between 1700 and 2100 CE.
%
% Function is vectorised.
%
% Input
% -----
% dtime = datetime in UTC time
% lat = local latitude (decimal degrees)
% lon = local longitude (decimal degrees)
% con = solar constant (leave empty, [], for 1361 W/m2)
%
% Output
% ------
% wm2 = irradiance (W/m2)
% sunlon = lambda, geocentric solar longitude (degrees)
% eot = equation of time (minutes)
% lappst = local apparent solar time (0-24)
% elev = solar elevation (degrees)
%
% B.C. Lougheed (bryan.lougheed@geo.uu.se)
% April 2023, Matlab 2019a

% --- manual function input for testing
% t1 = datetime('now','timezone','UTC');
% t2 = datetime('tomorrow','timezone','UTC');
% dtime = t1:minutes:t2;
% lat = 51.887511;
% lon = -8.4864234;
% con = 1361;

% check timezone input
if ~strcmpi('UTC',dtime.TimeZone)
	error('Please input dtime specified in UTC time, see Matlab datetime() timezone documentation for help')
end

if isempty(con)
	con = 1361;
end

% reshape dtime
dtime = reshape(dtime,numel(dtime),1);

% get eccentricity, obliquity and longitude of perihelion for the year
% Simon et al. (1994), Astron. Astrophys. 282, pp, 663-683
t = (year(dtime)-2000)/1000; % negative time is in the past
ecc = 0.0167086342-0.0004203654.*t-0.0000126734.*t.^2+1444e-10.*t.^3-2e-10.*t.^4+3e-10.*t.^5;
lpe = deg2rad( 102.93734808+(11612.35290.*t+53.27577.*t.^2+0.14095.*t.^3+0.1140.*t.^4+0.00478.*t.^5)/3600 );
obl = deg2rad( 23+(26/60)+(21.412/3600)-(468.0927.*t/3600)-(0.0152.*t.^2/3600)+(1.9989.*t.^3/3600)-(0.0051.*t.^4/3600)-(0.0025.*t.^5/3600) );

% calculate UTC datetime of spring equinox for the year
% Meeus (1998) Chapter 27
% Y = year(dtime)/1000;
% JDE0 = 1721139.29189 + 365242.13740*Y + 0.06134*Y^2 + 0.00111*Y^3 - 0.00071*Y^4;
% T = (JDE0 - 2451545.0) / 36525;
% W = 35999.373*T - 2.47;
% dlam = 1 + 0.0334*(cosd(W)) + 0.0007*(cosd(2*W));
% 
% S = 485*cosd(324.96+1934.136*T);
% USE NASA TABLE FOR NOW
d = readtable('NASA_AR5_equinox.txt');
nasayr = d.Var1;
nasadate = d.Var2;
nasatime = d.Var3;
nasastring = cell(size(nasayr));
for i = 1:numel(nasayr)
	nasastring{i} = [num2str(nasayr(i)),'/',nasadate{i},' ',nasatime{i}];
end
nasaeqdt = datetime(nasastring,'InputFormat','yyyy/MM/dd HH:mm','timezone','UTC');
[~, idx] = ismember(year(dtime), nasayr);
yeareqdt = nasaeqdt(idx);

% mean solar days elapsed since spring equinox
delapsed = datenum(dtime-yeareqdt);

% calculate solar longitude and eot
tottime = 365.2425;
if delapsed < 0
	delapsed = deplapsed + tottime;
end
[sunlon eot] = time2sunlon(delapsed,ecc,lpe,tottime,obl);
sunlon = deg2rad(sunlon);

% get local apparent solar time, including correcting for eot
[h m s] = hms(minutes(eot) + dtime); % eot = solartime - utc   -->   solartime = eot + utc
pmappst = h + m/60 + s/3600;  % apparent solar time at prime meridian (0-24)
lappst = pmappst + lon*4/60; % local apparent solar time (0-24). four minutes time offset per degree lon
lappst(lappst<0) = lappst(lappst<0) + 24;
lappst(lappst>24) = lappst(lappst>24) - 24;

% calculate irradiance
% local hour angle (-180 to +180 deg, midday = 0 deg)
hangle = (2*pi/24) * (lappst - 12); % radians
% declination of the sun, https://en.wikipedia.org/wiki/Position_of_the_Sun
dsun = asin(sin(obl) .* sin(sunlon));
% solar elevation, https://en.wikipedia.org/wiki/Solar_zenith_angle
elev = asin( sin(dsun).*sin(deg2rad(lat)) + cos(dsun).*cos(deg2rad(lat)).*cos(hangle) );

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

% for function output in degrees
sunlon = rad2deg(sunlon);
elev = rad2deg(elev);

end % end function
