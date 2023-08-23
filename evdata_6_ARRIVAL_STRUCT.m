% Script to cycle through data for each eventand build structures with all
% necessary data for each arrival of interest. 
% 
% Z. Eilon 2016
clear 


%% parameters
phases = {'P','S','SKS','PKP','PKS'};
compsaves = {{'Z'},{'T'},{'R'},{'Z'},{'R'}}; % needs to be cell in cell (may be >1 chan per phase)
resamprate = 10 ; % new, common sample rate
wind = [-200 200]; % seconds before and after arrival to save data for this arrival
overwrite = false;

% % project details
% dbname = 'EARdb';
% dbdir = '/Users/zeilon/Work/EastAfrica/EARdb/'; % include final slash

% project details
dbname = 'FRES_PILOT';
dbdir = '~/Dropbox/Work/FRES_PILOT/'; % include final slash

%% Preliminaries
startorid = 1151;
wd = pwd;
addpath('matguts')
cd(dbdir);
run([dbdir,dbname,'_startup.m']);
spd = 24*60*60; % seconds per day - to convert time in seconds to serial time

% point to right data dir 
if exist(datadir,'dir')~=7, datadir = regexprep(datadir,'data','data-1'); end
% event details
load([infodir,'/events'],'evinfo'); 
% station details
load([infodir,'/stations'],'stainfo');

%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%
for ip = 1:length(phases)
    phase = phases{ip};compsave = compsaves{ip};
    
for ie = startorid:evinfo.norids %1:evinfo.norids % loop on orids % got to 870 for S
    orid = evinfo.orids(ie);
    fprintf('\n Orid %.0f %s \n\n',orid,evinfo.evtimes_IRISstr{ie}(1:end-4))
    evdir = [num2str(orid,'%03d'),'_',char(evinfo.datestamp(ie,:)),'/'];
    datinfofile = [datadir,evdir,'_datinfo'];

    if ~exist([datinfofile,'.mat'],'file'), fprintf('No data at all for this event\n'), continue, end
    load(datinfofile)
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end
    nstas = length(datinfo);

    ofile = [data_eqar_dir,evdir,'_EQAR_',phase];
    if exist([ofile,'_Z.mat'],'file')
        if overwrite
            yn = 'y';
        else
            yn = input(sprintf('%s EQAR file exists, overwrite? [y/n] ',phase),'s');
        end
        if strcmp(yn,'y')
            delete([ofile,'_Z.mat'])
            delete([ofile,'_R.mat'])
            delete([ofile,'_T.mat'])
        else
            fprintf('ok, skipping\n')
            continue
        end
    end
    
    % RESULTS STRUCTURE
    eqar = struct('phase',phase,...
                  'sta',{datinfo.sta}','slat',[],'slon',[],'selev',[],'isobs',[],...
                  'gcarc',0,'seaz',0,'rayp',0,'pred_arrT',0,...
                  'tt',[],'datZ',[],'datR',[],'datT',[],'datH',[],'corZ',[],'samprate',resamprate);

    wlen = diff(wind)*resamprate; % window length in samples
    
    % LOOP ON STAS
    for is = 1:nstas
        fprintf('Station %-8s ',sprintf('%s...',datinfo(is).sta))
        % APPLY DATA QUALITY CONDITIONS
%         if datinfo(is).crap == 1; 
%             fprintf('crap, skipping\n'),continue 
%         end
        if ~datinfo(is).rmresp
            fprintf('must have response removed\n'), continue
        end
        
        % GET STATION + ARRIVAL INFO
        load([datadir,evdir,datinfo(is).sta,'.mat']); % load sta data for this evt
        
        if ~any(strcmp({data.phases.phase},phase))
            fprintf('No %s arrival at this sta\n',phase), continue
        end
        
        % station details
        statmp = struct2cell(data.station);
        [eqar(is).sta,eqar(is).slat,eqar(is).slon,eqar(is).selev] = deal(statmp{:});
        eqar(is).isobs = ~isempty(which_OBS(data.station.name));
        
        % station-event details
        eqar(is).gcarc = data.gcarc;
        eqar(is).seaz = data.seaz;

        % arrival details
        ip = find(strcmp({data.phases.phase},phase),1,'first');
        eqar(is).rayp = data.phases(ip).rayparameter;
        eqar(is).pred_arrT = data.phases(ip).artime;
        
        % check resampling
        if resamprate>data.samprate
            warning('resamprate would alias data - use smaller resamprate')
            continue
        end
            
        % GET DATA
        tt = [data.phases(ip).artime + wind(1)/spd:1./resamprate/spd:data.phases(ip).artime + wind(2)/spd];
        tt = tt(1:wlen);
        eqar(is).tt = tt;
        if any(strcmp(data.chans.component,'Z'))
            eqar(is).datZ = interp1(data.tt,data.dat(:,strcmp(data.chans.component,'Z')),tt);
        end
        if any(strcmp(data.chans.component,'H'))
            eqar(is).datH = interp1(data.tt,data.dat(:,strcmp(data.chans.component,'H')),tt);
        end
        
        if ~datinfo(is).NEZ
            fprintf(' not rotated\n'), continue
        else
            if data.chans.azimuth(strcmp(data.chans.component,'N')) ~=0, error('Bad instrument orientation\n'), end
            if ~any(strcmp(data.chans.component,'E')), continue; end
            if ~any(strcmp(data.chans.component,'N')), continue; end
            datN = interp1(data.tt,data.dat(:,strcmp(data.chans.component,'N')),tt);
            datE = interp1(data.tt,data.dat(:,strcmp(data.chans.component,'E')),tt);
            
            foraz = mod(data.seaz+180,360);
            eqar(is).datR =  datN*cosd(foraz) + datE*sind(foraz);
            eqar(is).datT = -datN*sind(foraz) + datE*cosd(foraz); % positive to right looking along foraz
        end
        
        fprintf('got data\n')    
    end % loop on stas
    
    if length([eqar.slat])<3
        fprintf('NOT ENOUGH STATIONS (%.0f) TO MAKE IT WORTH SAVING\n',length([eqar.slat]))
        continue
    end

    % SAVE
    for ic = 1:length(compsave)
        if ~exist([data_eqar_dir,evdir],'dir'), mkdir([data_eqar_dir,evdir]); end
        save([ofile,'_',compsave{ic},'.mat'],'eqar')
    end
%     save([ofile,'_Z.mat'],'eqar')
%     save([ofile,'_R.mat'],'eqar')
%     save([ofile,'_T.mat'],'eqar')
    
    copyfile([datinfofile,'.mat'],[data_eqar_dir,evdir,'_datinfo_',phase,'.mat'])

    
end% loop on orids
end% loop on phases
cd(wd);
