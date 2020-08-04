function  [tka ecc] = getlaskar2010(option, varargin)
% [tka ecc] = getlaskar2010(option)
%
% Open Laskar2010 eccentricity solution data files.
% Downloaded from http://vo.imcce.fr/insola/earth/online/earth/La2010/index.html
%
% Input
% =====
%
% option = 1	La2010a (solution a)
% option = 2	La2010b (solution b)
% option = 3	La2010c (solution c)
% option = 4	La2010d (solution d)
%
% Output
% ======
%
% tka = time in ka BP (negative years = future)
% ecc = eccentricity
%
% Optional
% ========
% 'slice',[tmin tmax]
% Specify desired time interval (in ka BP)
%
% Cite:
% =====
% Laskar, J., Fienga, A., Gastineau, M., Manche, H., 2011. 
% La2010: a new orbital solution for the long-term motion of the Earth. 
% A&A 532, A89. https://doi.org/10.1051/0004-6361/201116836
%
% -----------------------------------------
% B.C. Lougheed / May 4, 2020 / Matlab2019a

% stupidly long way that matlab processes optional vargin
p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = false;
p.FunctionName='getlaskar2010';
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
	d = load('La2010a_ecc3L.dat');
elseif option == 2
	d = load('La2010b_ecc3L.dat');
elseif option == 3
	d = load('La2010c_ecc3L.dat');
elseif option == 4
	d = load('La2010d_ecc3L.dat');
end

d(:,1) = (d(:,1)*-1)-(50/1000); % ka 1950
d = sortrows(d,1,'descend');

if min(slice) == max(slice)
	d = d(d(:,1) == slice(1),:);
else
	d = d(d(:,1) >= min(slice) & d(:,1) <= max(slice),:);
end

tka = d(:,1);
ecc = d(:,2);



