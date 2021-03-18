% recreate figure data from Huybers (2006) doi:10.1126/science.1125249

% solar sonstant to use
con = 1361;
% latitude to look at
lattarg = 65;
% melting threshold W/m2
thresh = 275; % what huybers uses


% get laskar orbital parameters
[tka ecc obl lpe]  = getlaskar2004(1, 'slice',[-0.5 2000.5]);

% N65 summer solstice W/m2
[n65sswm2, ~, ~, ~] = irrwm2(65, 90, con, ecc, obl, lpe);

% total summer J/m2
ndays = 365.2;
dayres = 0.1;
sdays = 0:dayres:ndays-dayres; % day 0 and day ndays are same day
n65Jm2thresh = NaN(size(tka));
n65meanirrthresh = NaN(size(tka));
secs = NaN(size(tka));
for i = 1:numel(tka)
	sunlons = sday2sunlon(sdays,ecc(i),lpe(i),ndays);
	[irrs, ~, ~, ~] = irrwm2(lattarg, sunlons, con, ecc(i), obl(i), lpe(i));
	n65meanirrthresh(i) = mean(irrs(irrs>=thresh));
	secs(i) = numel(sdays(irrs>=thresh)) * (24*60*60*dayres);
	n65Jm2thresh(i) = n65meanirrthresh(i) * secs(i); % W/m2 * s = J/m2
end
seasdays = secs / (60*60*24); % season day length (days with Wm2 greater than threshold)


% plot some stuff

xlims = [0 1000];

subplot(3,1,1)
plot(tka,n65meanirrthresh)
set(gca,'xdir','reverse')
xlim(xlims)
ylabel('Wm^{-2}')
xlabel('Age (ka)')

subplot(3,1,2)
plot(tka,seasdays)
set(gca,'xdir','reverse')
xlim(xlims)
ylabel('Days')
xlabel('Age (ka)')

subplot(3,1,3)
plot(tka,n65Jm2thresh/10^9)
set(gca,'xdir','reverse')
xlim(xlims)
ylabel('GJm^{-2}')
xlabel('Age (ka)')
