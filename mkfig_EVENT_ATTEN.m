% SCRIPT TO MAKE A FIGURE OF THE ATTENUATION (SEEN IN TIME AND FREQ DOMAIN)
% OF A SINGLE EVENT - uses plot_ATTEN_TandF_domain
% 
% Z. Eilon 2016

clear all
close all
% cd /Users/zeilon/Documents/MATLAB/CASC_atten
addpath('matguts')

% project details
dbname = 'EARdb';
dbdir = '~/Dropbox/Work/EARdb/'; % include final slash


%% parameters
phase = 'P';
component = 'Z'; %'Z', 'R', or 'T'
orid = 2193; %263; %269/263 for S %275;

snrmin = 4;
% latlims = [40 48]; % [44.2 48]
% lonlims = [-132 -110]; % [-132 -120]

ifsave   = false;

%% directories 
% % ANTELOPE DB DETAILS
% dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
% dbnam = 'cascBIGdb';
% % DATA DIRECTORY (top level)
% datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash


%% Preliminaries
cd(dbdir);
run([dbdir,dbname,'_startup.m']);
spd = 24*60*60; % seconds per day - to convert time in seconds to serial time
run([plotdir,'map_parameters.m'])

%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

% % GET EVENTS DATA
% db = dbopen([dbdir,dbnam],'r');
% dbor = dblookup_table(db,'origin');
% dbors = dbsubset(dbor,sprintf('orid == %.0f',orid));
% [elat,elon,edep,evtime,mag] = dbgetv(dbors,'lat','lon','depth','time','ms');
% dbclose(db);
% 
% fprintf('\nOrid %.0f, M%.1f on %s at [%.2f, %.2f], %.0f km deep\n\n',...
%     orid,mag,epoch2str(evtime,'%Y-%m-%d'),elat,elon,edep)
% evdir       = [num2str(orid,'%03d'),'_',epoch2str(evinfo.evtime(ie)),'%Y%m%d%H%M'),'/'];

%% More modern version
% event details
load([infodir,'/events'],'evinfo'); 
% station details
load([infodir,'/stations'],'stainfo');

% name files and directories
ie = find(evinfo.orids==orid);
evdir = [num2str(orid,'%03d'),'_',evinfo.datestamp(ie,:),'/'];
datinfofile = [datadir,evdir,'_datinfo_',phase];
arfile      = [datadir,evdir,'_EQAR_',phase,'_',component];

% check files exist
if ~exist([datinfofile,'.mat'],'file'), fprintf('No station mat files for this event\n'), end
if ~exist([arfile,'.mat'],'file'), fprintf('No arfile for this event + phase\n'), end

% load files
load(datinfofile) % loads datinfo stuctur/e
load(arfile)      % loads eqar structure

% options to skip
if isempty(datinfo), fprintf('No station mat files for this event\n'), end
if all([eqar.pred_arrT]==0), fprintf('No %s arrivals for this event\n',phase), return, end


%% ========================== START PLOTTING =============================

% QC
indgd = 1:size(eqar);
indgd(isnan([eqar(indgd).snr_wf])) = []; % kill nan traces
indgd([eqar(indgd).snr_wf]<snrmin) = []; % kill low snr traces
indgd([eqar(indgd).slat]<latlims(1)) = []; % kill low-lat traces
indgd([eqar(indgd).slat]>latlims(2)) = []; % kill high-lat traces
indgd([eqar(indgd).slon]<lonlims(1)) = []; % kill too-east traces
indgd([eqar(indgd).slon]>lonlims(2)) = []; % kill too-west traces
% indgd(cellfun('isempty',regexp({eqar(indgd).sta},'([J,F,M,G]*[0-9][0-9][A-C]$)'))) = []; % kill land
% indgd(~cellfun('isempty',regexp({eqar(indgd).sta},'FN'))) = [];% kill 'FN___' stas
% indgd(~cellfun('isempty',regexp({eqar(indgd).sta},'M0'))) = [];% kill 'FN___' stas
if length(indgd) < 2, fprintf('NO GOOD TRACES/ARRIVALS, skip...\n'), return, end

plot_ATTEN_TandF_domain_fAxP_EAR( eqar(indgd) )


%% spectral ratios
figure(3), set(gcf,'position',[400 400 800 450])
set(gca,'fontsize',14,'XTick',[0.05:0.05:0.25])
ylim([-2.5 1.5])
if ifsave
    save2pdf(3,sprintf('1evt_specratios_%s_%s_%.0f_%s',...
        phase,component,orid,epoch2str(evtime,'%Y-%m-%d')),'figs/1evt')
end

%% waveforms
figure(2), 
set(gca,'fontsize',14)
if ifsave
    save2pdf(2,sprintf('1evt_waveforms_%s_%s_%.0f_%s',...
        phase,component,orid,epoch2str(evtime,'%Y-%m-%d')),'figs/1evt')
end

%% mapview
figure(31), 
set(gcf,'position',[400 400 750 450])
axis([lonlims latlims]+[-1 0 -1.4 1.5])
set(gca,'fontsize',14)
title('$\Delta t^*$ recorded across JdF OBS stations','interpreter','latex','fontsize',18)

if ifsave
    save2pdf(31,sprintf('1evt_dtstar_map_%s_%s_%.0f_%s',...
        phase,component,orid,epoch2str(evtime,'%Y-%m-%d')),'figs/1evt')
end

%% section
figure(17), 
% title('Section of $\Delta t^*$ and $\delta T$ recorded across JdF ','interpreter','latex','fontsize',18)
if ifsave
    save2pdf(17,sprintf('1evt_section_%s_%s_%.0f_%s',...
        phase,component,orid,epoch2str(evtime,'%Y-%m-%d')),'figs/1evt')
end

