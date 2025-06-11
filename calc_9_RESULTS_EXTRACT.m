%% parse results into large nstas*nevts structures
clear all
% % % if running alone, establish these, otherwise use previous values
% % if ~exist('phase') || ~exist('component') 
% % phase = 'P';
% % component = 'Z';
% % end

% % project details
dbname = 'EARdb';
% % dbdir = '/Volumes/Lacie/Granite_EastAfrica/EARdb/'; % include final slash
dbdir = '/Users/zeilon/Dropbox/Work/EARdb/'; % include final slash

% project details
% dbname = 'FRES_PILOT';
% dbdir = '~/Dropbox/Work/FRES_PILOT/'; % include final slash

%% Preliminaries
wd = pwd;
addpath('matguts')
cd(dbdir);
run([dbdir,dbname,'_startup.m']);
spd = 24*3600;


%% Loop over all phases, for efficiency
phases = {'P','S','SKS','PKS','PP'};
components = {'Z','T','R','R','Z'};
% phases = {'PP'};
% components = {'Z'};

methods = {'dT','specR','fAxP'}; % choose from among {'dT','specR','comb','fAxP'};

for idata = 1:length(phases)
    phase = phases{idata};
    component = components{idata};

%% conditions
mag_min = 6.25; % skip events if below this mag
acor_min = 0.75; % skip events if xcor max below this acor
snr_min = 3; % skip result if data below this snr
event_dtime_skip = 60; % in seconds, skip event if less than this time from previous event (and used result for prev event)




% resdir = '/Volumes/LaCie/Granite_EastAfrica/EARdb/RESULTS_p2a_5/';

%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

% event details
load([infodir,'/events'],'evinfo'); 
% station details
load([infodir,'/stations'],'stainfo');

if any(strcmp(methods,'dT'))
    all_dT     = nan(stainfo.nstas,evinfo.norids);
end
if any(strcmp(methods,'specR'))
    all_dtstar_specR = nan(stainfo.nstas,evinfo.norids);
end
if any(strcmp(methods,'comb'))
    all_dT_comb = nan(stainfo.nstas,evinfo.norids);
    all_dtstar_comb = nan(stainfo.nstas,evinfo.norids);
    all_dtstar_comb_std = nan(stainfo.nstas,evinfo.norids);
end
if any(strcmp(methods,'fAxP'))
    all_dT_fAxP = nan(stainfo.nstas,evinfo.norids);
    all_dtstar_fAxP = nan(stainfo.nstas,evinfo.norids);
    all_dtstar_fAxP_std = nan(stainfo.nstas,evinfo.norids);
end

for ie = 1:evinfo.norids  % loop on orids
    if  evinfo.evmags(ie)<mag_min, continue, end

    % skip if a spurious catalog entry repeat of previous event
    if ie>1
    if abs(evinfo.evtimes(ie) - evinfo.evtimes(ie-1))*spd < event_dtime_skip % this is a repeat
        if any(isnan(all_dT(:,ie-1))), continue; end % and there is already data from first one
    end
    end
    fprintf('Orid %.0f\n',ie)

    orid = evinfo.orids(ie);
    % name files and directories
    evdir = [num2str(orid,'%03d'),'_',evinfo.datestamp(ie,:),'/'];
    datinfofile = [data_eqar_dir,evdir,'_datinfo_',phase];
    arfile      = [data_eqar_dir,evdir,'_EQAR_',phase,'_',component];
       % check files exist
    if ~exist([datinfofile,'.mat'],'file'), fprintf('No station mat files for this event\n');continue, end
    if ~exist([arfile,'.mat'],'file'), fprintf('No arfile for this event + phase\n');continue, end    
    % load info file
    ld = load(datinfofile,'datinfo'); datinfo = ld.datinfo; % loads datinfo stucture
    
    % options to skip
    if isempty(datinfo), fprintf('No station mat files for this event\n'); continue, end
    if ~isfield(datinfo,'xcor'), fprintf('   NEED TO XCOR\n',phase), continue, end

    % load data file
    ld = load(arfile); eqar = ld.eqar;      % loads eqar structure
    if all([eqar.pred_arrT]==0), fprintf('No %s arrivals for this event\n',phase), continue, end
    
%% parse dT data
if any(strcmp(methods,'dT'))
    if ~isfield(eqar,'dT'), fprintf('   NEED TO XCOR\n',phase), continue, end
      
    for is = 1:stainfo.nstas
        ise = find( strcmp({datinfo.sta},stainfo.stas(is)) & strcmp({datinfo.nwk},stainfo.nwk(is))  ); 
        if isempty(ise), continue; end
        if (eqar(ise).acor < acor_min) || isnan(eqar(ise).acor), continue; end
        if (isempty(eqar(ise).dT)) || isnan(eqar(ise).dT), continue, end
        all_dT(is,ie) = eqar(ise).dT;
    end
end
        
%% parse dtstar - mtm spec ratio data
if any(strcmp(methods,'specR'))
    if ~isfield(eqar,'dtstar_specR')
        if ~isfield(eqar,'dtstar')
            fprintf('   NEED TO CALC DTSTAR\n',phase)
            continue, 
        else
            warning('Seems to be calling specR dtstar just "dtstar"')
            try
                eqar.dtstar_specR = eqar.dtstar;
            catch
                eqar = dealto(eqar,'dtstar_specR',[eqar.dtstar]);
            end
        end
    end

    for is = 1:stainfo.nstas
        ise = find( strcmp({datinfo.sta},stainfo.stas(is)) & strcmp({datinfo.nwk},stainfo.nwk(is))  ); 
        if isempty(ise), continue; end
        if eqar(ise).snr_wf < snr_min, continue; end
        if isempty(eqar(ise).dtstar_specR), continue, end
        all_dtstar_specR(is,ie) = eqar(ise).dtstar_specR;
    end
end
    
%% parse dtstar_COMB data
if any(strcmp(methods,'comb'))
    nocomb = false;
    if ~isfield(eqar,'dtstar_comb'), fprintf('   NEED TO CALC DTSTAR w COMB\n',phase), nocomb = true; end
    if ~isfield(eqar,'dtstarstd_comb'), fprintf('   NEED TO REDO COMB FOR THIS EVT\n',phase), nocomb = true; end

    if ~nocomb
        for is = 1:stainfo.nstas
            ise = find( strcmp({datinfo.sta},stainfo.stas(is)) & strcmp({datinfo.nwk},stainfo.nwk(is))  ); 
            if isempty(ise), continue; end
            if eqar(ise).snr_wf < snr_min, continue; end
            if isempty(eqar(ise).dtstar_comb), continue, end
            if isnan(eqar(ise).dtstar_comb), continue, end
            all_dT_comb(is,ie) = eqar(ise).dT_comb;
            all_dtstar_comb(is,ie) = eqar(ise).dtstar_comb;
            all_dtstar_comb_std(is,ie) = eqar(ise).dtstarstd_comb;
        end
    end
end

%% parse dtstar_fAxP data
if any(strcmp(methods,'fAxP'))
    nofaxp = false;
    if ~isfield(eqar,'dtstar_fAxP'), fprintf('   NEED TO CALC DTSTAR w fAxP\n',phase), nofaxp = true; end
    if ~isfield(eqar,'dtstarstd_fAxP'), fprintf('   NEED TO REDO fAxP FOR THIS EVT\n',phase), nofaxp = true; end

    if ~nofaxp
        for is = 1:stainfo.nstas
            ise = find( strcmp({datinfo.sta},stainfo.stas(is)) & strcmp({datinfo.nwk},stainfo.nwk(is))  ); 
            if isempty(ise), continue; end
            if eqar(ise).snr_wf < snr_min, continue; end
            if isempty(eqar(ise).dtstar_fAxP), continue, end
            if isnan(eqar(ise).dtstar_fAxP), continue, end
            all_dT_fAxP(is,ie) = eqar(ise).dT_fAxP;
            all_dtstar_fAxP(is,ie) = eqar(ise).dtstar_fAxP;
            all_dtstar_fAxP_std(is,ie) = eqar(ise).dtstarstd_fAxP;
        end
    end
end
    
end % loop on orids

%% take off means of differential values
if any(strcmp(methods,'dT'))
    all_dT = all_dT - ones(stainfo.nstas,1)*nanmean(all_dT);
end
if any(strcmp(methods,'specR'))
    all_dtstar_specR = all_dtstar_specR - ones(stainfo.nstas,1)*nanmean(all_dtstar_specR);
end
if any(strcmp(methods,'comb'))
    all_dtstar_comb = all_dtstar_comb - ones(stainfo.nstas,1)*nanmean(all_dtstar_comb);
end
if any(strcmp(methods,'fAxP'))
    all_dtstar_fAxP = all_dtstar_fAxP - ones(stainfo.nstas,1)*nanmean(all_dtstar_fAxP);
end

%% save
fprintf('SAVING ALL RESULT FILES to %s\n',resdir)

% dT
if any(strcmp(methods,'dT'))
    save([resdir,'all_dT_',phase,'_',component],'all_dT') % xcorr differential TT
end

% dtstar-specR
if any(strcmp(methods,'specR'))
    save([resdir,'all_dtstarspecR_',phase,'_',component],'all_dtstar_specR') % mtm specR differential tstar
end

% dtstar-comb
if any(strcmp(methods,'comb'))
    % save([resdir,'all_dTcomb_',phase,'_',component],'all_dT_comb') % IGNORE THIS ONE - NOT GOOD! comb amp+phase spec differential TT
    save([resdir,'all_dtstarcomb_',phase,'_',component],'all_dtstar_comb') % comb amp+phase spec differential tstar
    save([resdir,'all_dtstarcombstd_',phase,'_',component],'all_dtstar_comb_std') % comb amp+phase spec differential tstar
end
% dtstar-fAxP

if any(strcmp(methods,'fAxP'))
    % save([resdir,'all_dTfAxP_',phase,'_',component],'all_dT_fAxP') % IGNORE THIS ONE - NOT GOOD! fAxP amp+phase spec differential TT
    save([resdir,'all_dtstarfAxP_',phase,'_',component],'all_dtstar_fAxP') % comb amp+phase spec differential tstar
    save([resdir,'all_dtstarfAxPstd_',phase,'_',component],'all_dtstar_fAxP_std') % comb amp+phase spec differential tstar
end
cd(wd);

end % loop on phases
    
    
    
    