function sday = sunlon2sday(sunlon,ecc,lpe,totdays)
% sday = sunlon2sday(sunlon,ecc,lpe,totdays)
%
% Given a particular eccentricity and precession, get solar day of tropical year 
% associated with a particular geocentric solar longitude, i.e. by accounting for 
% conservation of angular momentum during orbit (Keplar 2nd Law).
%
% Input
% =====
% sunlon  = Keplarian geocentric solar longitude in degrees (i.e., "v" relative to NH spring) 
%		    0 = NH Spring, 90 = NH Summer, 180 = NH Autumn, 270 = NH Winter
%		    Either 1 value (used as constant if other inputs are array), or an array of values.
% ecc     = eccentricity, array
% lpe     = longitude of perihelion (a.k.a omega-bar), array (radians)
% totdays = total solar days in the year, single value. Use empty, [], for 365.24.
%
% Output
% ======
% sday    = Solar day of tropical year (where Day 0 is northern spring equinox).
%		    array, same size as ecc and lpe
%
% B.C. Lougheed, June 2020
% Matlab 2019a
%
% See following for background, as well as comments in the script:
% Berger (1978). https://doi.org/10.1175/1520-0469(1978)035%3C2362:LTVODI%3E2.0.CO;2
% Meeus, J., (1998). Astronomical Algorithms, 2nd ed. Willmann-Bell, Inc., Richmond, Virginia. (specifically Chapter 30).
% Berger et al. (2010): doi: 10.1016/j.quascirev.2010.05.007

if isempty(totdays) == 1
	totdays = 365.24;
end

% change lpe from heliocentric to geocentric: Fig 1 and Appendix B in e.g. Berger et al. (2010). Also Berger (1978)
omega = lpe+pi; % add 180
omega(omega>=2*pi) = omega(omega>=2*pi) - 2*pi; % wrap to 360

% First, get solar day of northern spring equinox relative to perihelion
veq = 2*pi - omega; % v of spring equinox relative to perihelion
Eeq = 2 * atan( tan(veq/2) .* sqrt((1-ecc)./(1+ecc)) ); % Meeus (1998) page 195, solve for E
Meq = Eeq-ecc.*sin(Eeq); % Meeus page 195, solve for M (Keplar equation). M is the circular orbit equivalent of v
Meq(Meq<0) = pi + (pi-Meq(Meq<0)*-1); % inbound to perihelion
deq = Meq/(2*pi) .* totdays;

% Second, get solar day of target (x) sunlon relative to perihelion
vx = veq + deg2rad(sunlon);
vx(vx>2*pi) = vx(vx>2*pi) - 2*pi;
Ex = 2 * atan( tan(vx/2) .* sqrt((1-ecc)./(1+ecc)) ); % Meeus (1998) page 195, solve for E
Mx = Ex-ecc.*sin(Ex);% Meeus page 195, solve for M (Keplar equation). M is the circular orbit equivalent of v
Mx(Mx<0) = pi + (pi-Mx(Mx<0)*-1); % inbound to perihelion
dx = Mx/(2*pi) .* totdays;

% Finally, get solar day of target (x) sunlon relative to tropical year (NH spring equinox)
dx(dx<deq) = dx(dx<deq) + totdays; % for those in next orbital period relative to perihelion, keep in same orbital period relative to NH spring equinox
sday = dx - deq;

% eliminate rounding errors at zero
sday(sunlon==0) = 0;
sday(sunlon==360) = 0;


end % end main function

