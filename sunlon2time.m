function [time eot] = sunlon2time(sunlon,ecc,lpe,tottime,obl)
% [sday eot] = sunlon2time(sunlon,ecc,lpe,tottime,obl)
%
% Given a particular eccentricity and longitude of perihelion, get time of tropical year 
% associated with a particular geocentric solar longitude, i.e. by accounting for 
% conservation of angular momentum during orbit (Kepler 2nd Law).
%
% Input
% =====
% sunlon  = Keplerian geocentric solar longitude in degrees ('lambda', i.e. 'v' relative to NH spring equinox) 
%           0 = NH Spring, 90 = NH Summer, 180 = NH Autumn, 270 = NH Winter
%           Either 1 value (used as constant if other inputs are array), or a 1D array of values.
% ecc     = eccentricity, array
% lpe     = heliocentric longitude of perihelion (a.k.a omega-bar), array (in radians)
% tottime = total time in the year, single value, any constant time unit you want. Use empty, [], for 365.24.
% obl     = optional fifth input variable, for when calculating eot. Obliquity (in radians).
%
% Apologies that sunlon is in degrees and lpe in radians, it seemed like a
% good idea at the time.
%
% Output
% ======
% time    = Solar day of tropical year (where Day 0 is northern spring equinox).
%           Same dims as sunlon.
% eot     = equation of time (minutes). Requires fifth input variable obl.
%			Same dims as sunlon. Returns empty if obl not supplied.
%
% B.C. Lougheed, June 2020, Matlab 2019a
% Updated April 2023 to include eot.
%
% See following for background, as well as comments in the script:
% Berger (1978). https://doi.org/10.1175/1520-0469(1978)035%3C2362:LTVODI%3E2.0.CO;2
% Meeus, J., (1998). Astronomical Algorithms, 2nd ed. Willmann-Bell, Inc., Richmond, Virginia. (specifically Chapter 30).
% Berger et al. (2010): doi: 10.1016/j.quascirev.2010.05.007
% https://dr-phill-edwards.eu/Science/EOT.html (for equation of time)

if isempty(tottime) == 1
	tottime = 365.24;
end

% change lpe from heliocentric to geocentric: Fig 1 and Appendix B in e.g. Berger et al. (2010). Also Berger (1978)
omega = lpe+pi; % add 180
omega(omega>=2*pi) = omega(omega>=2*pi) - 2*pi; % wrap to 360

% First, get day of anchor day (dz) relative to perihelion
vz = 2*pi - omega; % v of spring equinox relative to perihelion
vz(vz>2*pi) = vz(vz>2*pi) - 2*pi;
Ez = 2 * atan( tan(vz/2) .* sqrt((1-ecc)./(1+ecc)) ); % Meeus (1998) page 195, solve for E
Mz = Ez-ecc.*sin(Ez); % Meeus page 195, solve for M (Kepler equation). M is the circular orbit equivalent of v
Mz(Mz<0) = pi + (pi-Mz(Mz<0)*-1); % inbound to perihelion
dz = Mz/(2*pi) .* tottime;

% Second, get day of target day (dx) relative to perihelion
vx = vz + deg2rad(sunlon);
vx(vx>2*pi) = vx(vx>2*pi) - 2*pi;
Ex = 2 * atan( tan(vx/2) .* sqrt((1-ecc)./(1+ecc)) ); % Meeus (1998) page 195, solve for E
Mx = Ex-ecc.*sin(Ex); % Meeus page 195, solve for M (Kepler equation). M is the circular orbit equivalent of v
Mx(Mx<0) = pi + (pi-Mx(Mx<0)*-1); % inbound to perihelion (probably not necessary)
dx = Mx/(2*pi) .* tottime;

% Finally, get day of target day (dx) relative to day of anchor day (dz)
dx(dx<dz) = dx(dx<dz) + tottime; % for dz in next orbital period relative to perihelion, keep in same orbital period relative to NH spring equinox
time = dx - dz;

% eliminate rounding errors at zero
time(sunlon==0) = 0;
time(sunlon==360) = 0;

% calculate eot
if nargin > 4
	sunlon = deg2rad(sunlon);
	% https://dr-phill-edwards.eu/Science/EOT.html (somebody who actually explains it clearly) 
	% eccentricity component
	dtecc = rad2deg(Mx-vx) * 4; % four minutes per degree
	% obliquity component
	alpha = atan2(sin(sunlon) * cos(obl), cos(sunlon));
	alpha(alpha<0) = alpha(alpha<0) + 2*pi;
	dtobl = rad2deg((sunlon-alpha)) * 4;
	eot = dtecc + dtobl;
else
	eot = [];	
end

end % end main function

