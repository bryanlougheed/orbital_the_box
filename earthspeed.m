function [v] = earthspeed(rx)
% [v] = earthspeed(rx)
%
% Input
% -----
% rx = distace from the Sun in AU. Single value or vector.
%
% Output
% ------
% v   = instantaneous velocity of the
%       Earth in its eccentric orbit
%       at moment when distance from Sun is rx.
%       Assuming semi-major axis of 1 AU.
%       Output is in m/s, same dims as rx
%
% B.C. Lougheed
% December 2020, Matlab 2020a
% See https://en.wikipedia.org/wiki/Orbital_speed
% and https://www.youtube.com/watch?v=X1XrtqubGPE

% I use very precise numbers here (from wikipedia, ha)
au = 1.495978707*10^11; % au in m
rx = rx * au;
mus = 1.32712440018*10^20; % grav constant of sun (ignoring earth... too tiny)

% rp = (1-ecc) * au; % perihelion distance
% ra = (1+ecc) * au; % aphelion distance
% a = (rp+ra)/2;

v = sqrt(mus .* ( 2./rx - 1/au));


