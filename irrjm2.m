function [intirr, ndays] = irrjm2(thresh, lat, con, dayres, ecc, obl, lpe, earthshape)
% [intirr, ndays] = irrjm2(thresh, lat, con, dayres, ecc, obl, lpe, earthshape)
%
% Calculate integrated irradiation (J/m2) at top of atmosphere for all days
% exceeding a certain threashold in mean irradiance (W/m2).
%
% Input
% =====
%
% thresh = Threshold value (W/m2). Single value, or vector of values.
% lat = Geocentric latitude (in deg. N, minus for S) on Earth. Single valule.
% con = Solar constant. Single numerical value or 1D array, W/m2. Leave empty, i.e. [], for 1361.
% dayres = day resolution for the integration, must be divisible into 365.2. Recommended to use 0.1.
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
% intirr = Calculated integrated irradiation at top of atmosphere for days exceeding thresh. 
%          J/m2. Array same dims as ecc, obl and lpe.
%
% ndays  = Number of days exceeding thresh. Same dims as intirr.
%
% -----------------------------------------
%
% B.C. Lougheed, August 2021, Matlab 2020a

if isempty(con)
	con = 1361;
end

if nargin < 8
	earthshape = 'sphere';
end

days = 0:dayres:365.2-dayres;

intirr = NaN(size(ecc));
ndays = NaN(size(ecc));
for i = 1:numel(ecc)
	sunlons = sday2sunlon(days,ecc(i),lpe(i),365.2);
	insos = irrwm2(lat, sunlons, con, ecc(i), obl(i), lpe(i), earthshape);
	ndays(i) = numel(days(insos>=thresh)) * dayres;
	intirr(i) = mean(insos(insos>=thresh)) * (ndays(i)*24*60*60); % W/m2 to J/m2
end

end