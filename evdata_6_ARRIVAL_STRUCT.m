% Script to cycle through data for each eventand build structures with all
% necessary data for each arrival of interest. 
% 
% Z. Eilon 2016
clear 


%% parameters
phases = {'P','PP','S','SKS','PKP','PKS'};
compsaves = {{'Z'},{'Z'},{'T'},{'R'},{'Z'},{'R'}}; % needs to be cell in cell (may be >1 chan per phase)
phases = {'PP'};
compsaves = {{'Z'}}; % needs to be cell in cell (may be >1 chan per phase)
resamprate = 10 ; % new, common sample rate
wind = [-200 200]; % seconds before and after arrival to save data for this arrival

overwrite = true; % overwrite all
overwrite_blank = true; % overwrite if no measurements for a given eqar.
overwrite_Dsta_above = 3; % overwrite if the number of new stations is this value or higher. Leave 0 or false to ignore

% % project details
dbname = 'EARdb';
dbdir = '/Users/zeilon/Dropbox/Work/EARdb/'; % include final slash

% project details
% dbname = 'FRES_PILOT';
% dbdir = '~/Dropbox/Work/FRES_PILOT/'; % include final slash

%% Preliminaries
startorid = 1;
% time bounds for events to get - ignore before startdate, or after enddate
startdate = '0204-03-01'; % format 'YYYY-MM-DD' 
enddate   = '2025-01-01'; % format 'YYYY-MM-DD'

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
    
for ie = startorid:evinfo.norids %1:evinfo.norids % loop on orids % got to 870 for S
    orid = evinfo.orids(ie);

    % ignore outside date bounds
    if evinfo.evtimes(ie) < datenum(startdate) || evinfo.evtimes(ie) > datenum(enddate), continue; end

    fprintf('\n Orid %.0f %s \n\n',orid,evinfo.evtimes_IRISstr{ie}(1:end-4))
    evdir = [num2str(orid,'%03d'),'_',char(evinfo.datestamp(ie,:)),'/'];
    datinfofile = [datadir,evdir,'_datinfo'];

    if ~exist([datinfofile,'.mat'],'file'), fprintf('No data at all for this event\n'), continue, end
    load(datinfofile)
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end
    nstas = length(datinfo);

for ipp = 1:length(phases)
    phase = phases{ipp}; compsave = compsaves{ipp};
    fprintf('>> phase %s\n',phase);

    ofile = [data_eqar_dir,evdir,'_EQAR_',phase];
    if exist([ofile,'_Z.mat'],'file') || exist([ofile,'_R.mat'],'file') || exist([ofile,'_T.mat'],'file')
        clear('yn')
        if overwrite
            yn = 'y';
        else
            fprintf('%s EQAR file exists\n',phase)
            datinfo_phase = load([data_eqar_dir,evdir,'_datinfo_',phase,'.mat']);
            if ~isfield(datinfo_phase.datinfo,'xcor')
                fprintf('   BUT seems to have no dT (xcor) data\n')
                if overwrite_blank
                    yn = 'y';
                else
                    yn = input('Overwrite EQAR?? [y/n] ','s');
                end
            else
                fprintf('!! Already has **%.0f** xcor data\n',sum([datinfo_phase.datinfo.xcor]))
                fprintf('!! BUT old file has %.0f stations, whereas now have %.0f stations\n',...
                    length({datinfo_phase.datinfo.sta}),length({datinfo.sta}))
                fprintf('!! WARNING: OVERWRITE WILL RESET MEASUREMENTS FOR THIS EVENT/PHASE\n')
                
                if overwrite_Dsta_above && ...
                        length({datinfo.sta}) - length({datinfo_phase.datinfo.sta}) >= overwrite_Dsta_above
                    yn = 'y';
                else
                    yn = input('Overwrite EQAR?? [y/n] ','s');
                end
                % since already some data, save it for safekeeping
                if any(strcmp(yn,{'y','yes','Y','YES'}))
                    ydaystr = datestr(now-1,'YYYYmmDD');
                    for icomp = ['ZRT']
                    if exist([ofile,'_',icomp,'.mat'],'file')
                        copyfile([ofile,'_',icomp,'.mat'],[ofile,'_',icomp,'_',ydaystr,'.mat']); 
                    end
                    end
                end
            end
            
        end
        if any(strcmp(yn,{'y','yes','Y','YES'}))
            fprintf('--OVERWRITING--\n')
            if exist([ofile,'_Z.mat'],'file'), delete([ofile,'_Z.mat']); end
            if exist([ofile,'_R.mat'],'file'), delete([ofile,'_R.mat']); end
            if exist([ofile,'_T.mat'],'file'), delete([ofile,'_T.mat']); end
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
        
        % GET phase timing from data structure. 
        %  if not there, use TauP. 
        %  if no arrival, move on.
        if ~any(strcmp({data.phases.phase},phase))
            % phase and distance details
            TT = tauptime('event',[evinfo.elats(ie),evinfo.elons(ie)],'depth',evinfo.edeps(ie),...
                            'station',[data.station.slat,data.station.slon],'phases',phase);
            if isempty(TT) % carry on if this phase not at this station 
                fprintf('No %s arrival at this sta\n',phase), continue
            else
                % only first instance
                TT = TT(1);
                TT.artime = TT.time/spd + evinfo.evtimes(ie);     
                % insert into data.phases
                ipnew = length(data.phases)+1;
                TTflds = fieldnames(TT);
                for iTTfld = 1:length(TTflds)
                    data.phases(ipnew).(TTflds{iTTfld}) = TT.(TTflds{iTTfld});
                end
            end
        end
        
        % station details
        statmp = struct2cell(data.station);
        [eqar(is).sta,eqar(is).slat,eqar(is).slon,eqar(is).selev] = deal(statmp{:});
        eqar(is).isobs = ~isempty(which_OBS(data.station.name));
        
        % station-event details
        eqar(is).gcarc = data.gcarc;
        eqar(is).seaz = data.seaz;

        % arrival details
        ipi = find(strcmp({data.phases.phase},phase),1,'first');
        eqar(is).rayp = data.phases(ipi).rayparameter;
        eqar(is).pred_arrT = data.phases(ipi).artime;
        
%         % check resampling
%         if resamprate < data.samprate % downsampling a bad idea; need to put on anti-alias filter
%             warning('resamprate would alias data - use smaller resamprate')
%             continue
%         elseif resamprate > data.samprate
%             % CAREFUL - YOU ARE RESETTING THE DATA HERE. BUT NOT SAVING!
%             data.dat = resample(data.dat,resamprate,data.samprate); 
%             data.tt = data.tt(1) + [0:(size(data.dat,1)-1)]'./resamprate./spd;
%         end

% better resampling - interpolate if upsamp, anti-alias filter if downsamp.
        % CAREFUL - YOU ARE RESETTING THE DATA HERE. BUT NOT SAVING!
        data.dat = downsamp(data.dat,data.samprate,resamprate);
        data.tt = downsamp(data.tt,data.samprate,resamprate);
            
        % GET DATA
        tt = [data.phases(ipi).artime + wind(1)/spd:1./resamprate/spd:data.phases(ipi).artime + wind(2)/spd];
        tt = tt(1:wlen); 
        tt = tt(:); % for the love of god make this a column.
        eqar(is).tt = tt;
        if any(strcmp(data.chans.component,'Z'))
            eqar(is).datZ = interp1(data.tt,data.dat(:,strcmp(data.chans.component,'Z')),tt);
        end
        if any(strcmp(data.chans.component,'H'))
            eqar(is).datH = interp1(data.tt,data.dat(:,strcmp(data.chans.component,'H')),tt);
        end
        
        if ~datinfo(is).NEZ
            fprintf('not rotated (got vertical only!)\n'), continue
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
    
    if length([eqar.slat])<4
        fprintf('NOT ENOUGH STATIONS (%.0f) TO MAKE IT WORTH SAVING\n',length([eqar.slat]))
        continue
    end


    % SAVE
    fprintf('--SAVING--\n\n')    
    for ic = 1:length(compsave)
        if ~exist([data_eqar_dir,evdir],'dir'), mkdir([data_eqar_dir,evdir]); end
        save([ofile,'_',compsave{ic},'.mat'],'eqar')
    end
    
    copyfile([datinfofile,'.mat'],[data_eqar_dir,evdir,'_datinfo_',phase,'.mat'])

end% loop on phases (ipp)
    
end% loop on orids (ie)
cd(wd);
