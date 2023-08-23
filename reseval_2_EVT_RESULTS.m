function reseval_2_EVT_RESULTS(dbname,dbdir,phases,components)
% reseval_2_EVT_RESULTS(dbname,dbdir,phases,components)

arguments
    dbname
    dbdir = ['/Users/zeilon/Dropbox/Work/',dbname,'/'];
    phases = {'P','S','SKS'};
    components = {'Z','T','R'};
end


%% Evaluate all events and see what has results where

%% conditions
% mag_min = 6.25; % skip events if below this mag
% acor_min = 0.75; % skip events if xcor max below this acor
% snr_min = 10; % skip result if data below this snr

% project details
% dbname = 'EARdb';
% dbdir = '/Users/zeilon/Dropbox/Work/EARdb/'; % include final slash

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
evinfo.Nres = zeros(evinfo.norids,1);
for ip = 1:length(phases)
    phase = phases{ip};
    component = components{ip};

% load results
load([resdir,'all_dT_',phase,'_',component,'.mat']) % xcorr differential TT
load([resdir,'all_dtstarspecR_',phase,'_',component,'.mat']) % mtm specR differential tstar
% load([resdir,'all_dTcomb_',phase,'_',component,'.mat']) % IGNORE THIS ONE - NOT GOOD! combb amp+phase spec differential TT
load([resdir,'all_dtstarcomb_',phase,'_',component,'.mat']) % comb amp+phase spec differential tstar

evinfo.(['Nres_',phase,component,'_NTT']) = sum(~isnan(all_dT),1)';
evinfo.(['Nres_',phase,component,'_NtstR']) = sum(~isnan(all_dtstar_specR),1)';
evinfo.(['Nres_',phase,component,'_NtstC']) = sum(~isnan(all_dtstar_comb),1)';

evinfo.Nres = evinfo.Nres + sum(~isnan([all_dT;all_dtstar_specR;all_dtstar_comb]),1)';

% find station locations
for ifld = {'all_dT','all_dtstar_specR','all_dtstar_comb'}
for ie = 1:evinfo.norids
    eval(sprintf('inds = ~isnan(%s(:,ie));',ifld{:}))
    if ~any(inds), continue; end
    evinfo.(['Centroid_',phase,component,'_',fliplr(strtok(fliplr(ifld{:}),'_'))])(ie,:) ...
        = mean([stainfo.slons(inds),stainfo.slats(inds)]);
end
end


end

%% summary table
%header
fprintf(' Orid Date     Mag  LAT    LON    Nres')
for ii=1:3
for ip = 1:length(phases)
fprintf('  %s_%s',dtype{ii},phases{ip});
end
end
fprintf('\n\n');

for ie = 1:evinfo.norids
    if evinfo.Nres(ie)==0; continue; end
    
    fprintf('%4.0f %4s  %3.1f %8.4f %9.4f  %3u',evinfo.orids(ie),evinfo.evtimes_IRISstr{ie},evinfo.evmags(ie),evinfo.elats(ie),evinfo.elons(ie),evinfo.Nres(ie))
    for ii=1:3
    for ip = 1:length(phases)
        fprintf('   %3u',evinfo.(['Nres_',phases{ip},components{ip},'_N',dtype{ii}])(ie));
    end
    end
    fprintf('\n');
    
end
