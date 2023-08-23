%% Evaluate all stations and see what has results where
phases = {'P','S','SKS'};
components = {'Z','T','R'};

%% conditions
% mag_min = 6.25; % skip events if below this mag
% acor_min = 0.75; % skip events if xcor max below this acor
% snr_min = 10; % skip result if data below this snr

% project details
dbname = 'EARdb';
dbdir = '/Users/zeilon/Dropbox/Work/EARdb/'; % include final slash

%% Preliminaries
wd = pwd;
addpath('matguts')
cd(dbdir);
run([dbdir,dbname,'_startup.m']);

%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

% event details
load([infodir,'/events'],'evinfo'); 
% station details
load([infodir,'/stations'],'stainfo');

dtype = {'TT','tstR','tstC'};

%% grab and total all results
stainfo.Nres = zeros(stainfo.nstas,1);
for ip = 1:length(phases)
    phase = phases{ip};
    component = components{ip};

% load results
load([resdir,'all_dT_',phase,'_',component,'.mat']) % xcorr differential TT
load([resdir,'all_dtstarspecR_',phase,'_',component,'.mat']) % mtm specR differential tstar
load([resdir,'all_dTcomb_',phase,'_',component,'.mat']) % IGNORE THIS ONE - NOT GOOD! combb amp+phase spec differential TT
load([resdir,'all_dtstarcomb_',phase,'_',component,'.mat']) % comb amp+phase spec differential tstar

stainfo.(['Nres_',phase,component,'_NTT']) = sum(~isnan(all_dT),2);
stainfo.(['Nres_',phase,component,'_NtstR']) = sum(~isnan(all_dtstar_specR),2);
stainfo.(['Nres_',phase,component,'_NtstC']) = sum(~isnan(all_dtstar_comb),2);

stainfo.Nres = stainfo.Nres + sum(~isnan([all_dT,all_dtstar_specR,all_dtstar_comb]),2);
end

%% summary table
%header
fprintf(' STA   NWK    LAT    LON    Nres')
for ii=1:3
for ip = 1:length(phases)
fprintf('  %s_%s',dtype{ii},phases{ip});
end
end
fprintf('\n\n');

for is = 1:stainfo.nstas
    if stainfo.Nres(is)==0; continue; end
    
    fprintf('%4s  %2s %8.4f %9.4f  %3u',stainfo.stas{is},stainfo.nwk{is},stainfo.slats(is),stainfo.slons(is),stainfo.Nres(is))
    for ii=1:3
    for ip = 1:length(phases)
        fprintf('   %3u',stainfo.(['Nres_',phases{ip},components{ip},'_N',dtype{ii}])(is));
    end
    end
    fprintf('\n');
    
end
