% script to plot the comb of filters used for the calculation of spectra
% 
% Z. Eilon 2016

clear all
close all

Tmin = 1;
Tmax = 20;
Nwds = 25;
Tw_opt = 'scale';
npol = 4;
dt = 1./40;

ifsave = 0;
odir = 'figs';

%% prepare filter + cleaning parms
% Make set of period windows for bandpass filter
Tmids = logspace(log10(Tmin),log10(Tmax),Nwds)';
if strcmp(Tw_opt,'scale')
    Twdhs = 0.5*diff(logspace(log10(Tmin/2),log10(2*Tmax),Nwds+1)');
else
    Twdhs = 0.5*Tw_opt(:).*ones(size(Tmids));
end

% Twdhs = 0.5*Tmids;
fmids = 1./Tmids;

flos = 1./(Tmids + Twdhs);
fhis = 1./(Tmids - Twdhs);

%% plot it!
plot_filter_comb( flos,fhis,dt,npol,'F' )

set(gcf,'pos',[360   212   991   493])
% set(gca,'xscale','log','xlim',[1./40, 1.5])

%% save it
if ifsave
    save2pdf(33,'filter_comb',odir)
end