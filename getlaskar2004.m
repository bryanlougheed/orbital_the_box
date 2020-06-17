function  [tka ecc obl lpe pre] = getlaskar2004(option, varargin)
% [tka ecc obl lpe pre]  = getlaskar2004(option)
%
% Open Laskar2010 solution data files.
% Downloaded from http://vo.imcce.fr/insola/earth/online/earth/La2010/index.html
%
% Input
% =====
%
% option = 1	51 Ma to 0 Ma
% option = 2	0 Ma to 21 Ma in the future
% option = 3	101 Ma to 0 Ma
% option = 4	249 Ma to 0 Ma
% option = 5	51 Ma to 21 Ma in the future (concatenate options 1 & 2)
%
% Output
% ======
%
% tka = time in ka BP 1950 (negative years = future)
% ecc = eccentricity
% obl = obliquity (radians)
% lpe = longitude of perihelion from moving equinox (radians)
% pre = precession index, calculated as ecc*sin(per)
%
% Optional
% ========
% 'slice',[tmin tmax]
% Specify desired time interval (in ka BP)
%
% Cite:
% =====
% Laskar, J., Robutel, P., Joutel, F., Gastineau, M., Correia, A.C.M., Levrard, B., 2004.
% A long-term numerical solution for the insolation quantities of the Earth. 
% A&A 428, 261–285. https://doi.org/10.1051/0004-6361:20041335
%
% -----------------------------------------
% B.C. Lougheed / May 4, 2020 / Matlab2019a

% stupidly long way that matlab processes optional vargin
p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = false;
p.FunctionName='getlaskar2004';
defaultslice = [-inf inf];
if exist('OCTAVE_VERSION', 'builtin') ~= 0
	addParamValue(p,'slice',defaultslice,@isnumeric);
else
	if datenum(version('-date'))>datenum('May 19, 2013')
		addParameter(p,'slice',defaultslice,@isnumeric);
	else
		addParamValue(p,'slice',defaultslice,@isnumeric);
	end
end
parse(p,varargin{:});
slice = p.Results.slice;


% load data
if option == 1
	d = load('INSOLN.LA2004.BTL.ASC');
elseif option == 2
	d = load('INSOLP.LA2004.BTL.ASC');
elseif option == 3
	d = load('INSOLN.LA2004.BTL.100.ASC');
elseif option == 4
	d = load('INSOLN.LA2004.BTL.250.ASC');
elseif option == 5
	d1 = load('INSOLN.LA2004.BTL.ASC');
	d2 = load('INSOLP.LA2004.BTL.ASC');
	d = [d1(2:end,:); d2]; % both have a 0 ka step. identical 0 ka data in both.
end

d(:,1) = (d(:,1)*-1)-(50/1000); % ka 1950

if min(slice) == max(slice)
	d = d(d(:,1) == slice(1),:);
else
	d = d(d(:,1) >= min(slice) & d(:,1) <= max(slice),:);
end

tka = d(:,1);
ecc = d(:,2);
obl = d(:,3);
lpe = d(:,4);
pre = ecc.*sin(lpe);



