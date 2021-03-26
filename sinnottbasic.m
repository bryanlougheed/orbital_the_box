function E = sinnottbasic(M,ecc)
% E = sinnottbasic(M,ecc)
%
% Input
% -----
% M = mean anomaly (radians)
% ecc = eccentricity of the ellipse
%
% Output
% ------
% E = eccentric anomaly (radians)
%
% Roger Sinnott (1985) BASIC script, as suggested by Meeus (1998).
% Solves Kepler equation (e.g. Eq. 30.5 on Page 195 of Meeus, 1998) for E.
% First port to Matlab by Tiho Kostadinov (Kostadinov and Gilb, 2014). 
% Vectorised using logical indexing to process arrays by Bryan Lougheed in June 2020. 
%
% R.W. Sinnott (1985), "A computer assault on Kepler's equation", Sky and Telescope, vol. 70, page 159.
% Meeus, J., (1998). Chapter 30 in Astronomical Algorithms, 2nd ed. Willmann-Bell, Inc., Richmond, Virginia.
% Kostadinov and Gilb, (2014): doi:10.5194/gmd-7-1051-2014
F = sign(M);
M = abs(M)./(2*pi);
M = (M-floor(M))*(2*pi).*F;
M(M<0) = M(M<0)+(2*pi); % put in same relative orbit
F = ones(size(M));
F(M>pi) = -1;
M(M>pi) = (2*pi)-M(M>pi); % inbound
Eo = (pi/2)*ones(size(ecc));
D = (pi/4)*ones(size(ecc));
for j = 1:54 % Sinnot says number of iterations is 3.30 * significant figures of system. Matlab double has 16 digit precision, so 16*3.30=53. Let's use 54
	M1 = Eo-ecc.*sin(Eo);
	Eo = Eo + D.*sign(M-M1);
	D = D./2;
end
E = Eo.*F;


end

