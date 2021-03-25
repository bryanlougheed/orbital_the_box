function sunlon = sday2sunlon(sday,ecc,lpe,totdays)
% sunlon = sday2sunlon(sday,ecc,lpe,totdays)
%
% Given a particular eccentricity and longitude of perihelion, get geocentric solar longitude 
% associated with a particular solar day of tropical year i.e. by accounting for 
% conservation of angular momentum during orbit (Kepler 2nd Law).
%
% Input
% =====
% sday    = Solar day of tropical year (where Day 0 is northern spring equinox).
%	    One value or 1D array.
% ecc     = eccentricity, array
% lpe     = heliocentric longitude of perihelion (a.k.a. omega-bar), array (radians)
% totdays = total solar days in the year, single value. Use empty, [], for 365.24.
%
% Output
% ======
% sunlon  = Keplerian geocentric solar longitude in degrees (lambda, i.e. 'v' relative to NH spring equinox) 
%           0 = NH Spring, 90 = NH Summer, 180 = NH Autumn, 270 = NH Winter
%           Either 1 value (used as constant if other inputs are array), or a 1D array of values.
%           Same dims as sday.
%
% B.C. Lougheed, June 2020
% Matlab 2019a
%
% See following for background, as well as comments in the script:
% Berger (1978). https://doi.org/10.1175/1520-0469(1978)035%3C2362:LTVODI%3E2.0.CO;2
% R.W. Sinnott (1985), "A computer assault on Kepler's equation." Sky and Telescope, vol. 70, page 159.
% Meeus, J., (1998). Astronomical Algorithms, 2nd ed. Willmann-Bell, Inc., Richmond, Virginia. (specifically Chapter 30).
% Berger et al. (2010): doi: 10.1016/j.quascirev.2010.05.007
% Kostadinov and Gilb, (2014): doi: 10.5194/gmd-7-1051-2014

if isempty(totdays) == 1
	totdays = 365.24;
end

% change lpe from heliocentric to geocentric: Fig 1 and Appendix B in e.g. Berger et al. (2010) and also Berger (1978).
omega = lpe+pi; % add 180
omega(omega>=2*pi) = omega(omega>=2*pi) - 2*pi; % wrap to 360

% First, get variables of NH spring equinox relative to perihelion
veq = 2*pi - omega; % NH spring v relative to perihelion
Eeq = 2 * atan( tan(veq/2) .* sqrt((1-ecc)./(1+ecc)) ); % Meeus (1998) page 195, solve for E
Meq = Eeq-ecc.*sin(Eeq); % Meeus page 195, solve for M (Kepler equation). M is the circular orbit equivalent of v
Meq(Meq<0) = pi + (pi-Meq(Meq<0)*-1); % inbound to perihelion
deq = Meq/(2*pi) .* totdays; % day of spring equinox relative to perihelion

% Second, get v of target (x) v relative to perihelion
dx = deq + sday;
Mx = (dx./totdays) * 2*pi;
Ex = sinnottbasic(Mx,ecc); % solve Kepler equation for E, Sinnott (1985)
vx = 2 * atan( tan(Ex/2) .* sqrt((1+ecc)./(1-ecc)) ); % Meeus (1998) page 195, solve for v
vx(vx<0) = pi + (pi-vx(vx<0)*-1); % inbound to perihelion

% target day's v relative to NH spring equinox v
vx(vx<veq) = vx(vx<veq) + 2*pi;
sunlon = rad2deg(vx - veq);

% eliminate rounding errors at 0
sunlon(sday==0) = 0;
sunlon(sday==totdays) = 0;






end 

