% Script to build body wave dataset
% 
% Uses antelope origin table and site(chan) tables, looping through and
% then using the IRISrequest tools to get and store the data in a directory
% tree where each event has a folder containing .mat files that are the
% data for each station.
% 
% N.B. this could be significantly sped up by using a different request
% format/algorithm, but this one gives maximal control/certainty about what
% data you are getting.
% 
% N.B.2 Make sure to have allocated sufficient Java Heap Memory for this,
% or the request will fail sometimes. At least 1000 MB.
%  (go to MATLAB Preferences > MATLAB > General > Java Heap Memory)
% 
% Z. Eilon
clear 

% % project details
% dbname = 'EARdb';
% dbdir = '/Users/zeilon/Work/EastAfrica/EARdb/'; % include final slash

% project details
dbname = 'FRES_PILOT';
dbdir = '~/Dropbox/Work/FRES_PILOT/'; % include final slash

datawind = [-100 1700]; % time window in seconds after event to [start end]

phases = 'P,S,PKP,PKS,SKS';

overwrite = false;

resamprate = 10; % leave empty or zero to use existing samprate

getnoise = false;
% OBSnoiseprewind = [-43200 0]; % time window in seconds after event to [start end] <== 12 hours in advance


%% Preliminaries
wd = pwd;
addpath('matguts')
cd(dbdir);
run([dbdir,dbname,'_startup.m']);
javaaddpath('IRIS-WS-2.0.15.jar')
spd = 24*60*60; % seconds per day - to convert time in seconds to serial time

% point to right data dir 
if exist(datadir,'dir')~=7, datadir = regexprep(datadir,'data','data-1'); end
% event details
load([infodir,'/events'],'evinfo'); 
% station details
load([infodir,'/stations'],'stainfo');

for ie = 1979:evinfo.norids % got to 3691 % 160-335 is Y6
    % sort out event stuff
    orid = evinfo.orids(ie);
    elat = evinfo.elats(ie); elon = evinfo.elons(ie); edep = evinfo.edeps(ie); 
    evtime = evinfo.evtimes(ie);
    evdir = [num2str(orid,'%03d'),'_',evinfo.datestamp(ie,:)];
    
    % make event directory
    if exist([datadir,evdir],'dir')~=7, mkdir(datadir,evdir); end
    
    % datinfo
    datinfofile = [datadir,'/',evdir,'/_datinfo'];
    if exist([datinfofile,'.mat'],'file')~=2 || overwrite==true
        datinfo = struct('sta',[],'nwk',[],'chans',[],'NEZ',false,'rmresp',false,'rmtilt',false,'rmcomp',false,'spectra',false);
    else
        load(datinfofile);
    end
        
    
    fprintf('REQUESTING DATA FOR EVENT %.0f (%s)\n',orid,evdir)
    for is = 1:stainfo.nstas % 1:nstas
        sta = stainfo.stas{is};
        datafile = [datadir,'/',evdir,'/',stainfo.stas{is}];
        
        if overwrite==false && exist([datafile,'.mat'],'file')==2, fprintf('   data already exists for %s %s\n',stainfo.nwk{is},sta); continue; end
        
%         % this is where we'll pull out the channel and station on/off info
%         db = dbopen([dbdir,dbnam],'r');
%         dbsch = dblookup_table(db,'sitechan');
%         dbschs = dbsubset(dbsch,sprintf('sta == "%s"',stas{is}));
%         nchans = dbnrecs(dbschs);
%         if nchans==0, dbclose(db); continue; end
%         [chans,ondates,offdates,azimuths,dips] = dbgetv(dbschs,'chan','ondate','offdate','hang','vang');
%         dbclose(db);
%         if ~iscell(chans), chans = {chans}; end;
%         
        % calc. data window
        if stainfo.selevs(is)<0 % if OBS, if the option is selected, grab big noise window too!
            if getnoise
                STARTtime = evtime + datawind(1)/spd + OBSnoiseprewind(1)/spd - 1/spd; % inc buffer
                ENDtime   = evtime + datawind(2)/spd + OBSnoiseprewind(2)/spd + 1/spd; % inc buffer
            else
                STARTtime = evtime + datawind(1)/spd - 1/spd; % inc buffer
                ENDtime   = evtime + datawind(2)/spd + 1/spd; % inc buffer
            end
        else % if LAND
        STARTtime = evtime + datawind(1)/spd - 1/spd; % inc buffer
        ENDtime   = evtime + datawind(2)/spd + 1/spd; % inc buffer
        end
        
        % check this station was alive for this event
        if stainfo.ondate(is)  > STARTtime , continue, end % continue if stat turned on after datwind start
        if stainfo.offdate(is) < ENDtime, continue, end % continue if stat turned off before datwind end

        STARTtime_str = datestr(STARTtime,'yyyy-mm-dd HH:MM:SS.FFF');
        ENDtime_str   = datestr(ENDtime,'yyyy-mm-dd HH:MM:SS.FFF');
        
        % make string list of chans to request
        chans = stainfo.chans(is,1:stainfo.nchans(is));
        chreq = []; for ic = 1:length(chans), chreq = [chreq,chans{ic},',']; end; chreq = chreq(1:end-1); %#ok<AGROW>
        
%         fprintf('   request station %.0f %s... ',is,stas{is})
        fprintf('   request station %.0f %s %s %s (%s - %s) ',is,stainfo.nwk{is},stainfo.stas{is},chreq,STARTtime_str,ENDtime_str)
        
        %% ======== GET THE DATA =======
        clear trace; % to be really, really sure
%         trace=irisFetch.Traces(nwk{is},sta_use,'*',chreq,STARTtime_str,ENDtime_str,'VERBOSE');
        trace=irisFetch.Traces(stainfo.nwk{is},stainfo.stas{is},'*',chreq,STARTtime_str,ENDtime_str);
        
        if isempty(trace), fprintf('NO DATA\n'); continue; end
        [ trace ] = fixtrace( trace );
        fprintf('got chans '); for ic = 1:length(trace), fprintf('%s, ',trace(ic).channel); end
        % if B and H channels - only keep H (might need higher samprate?)
        if strcmp([trace.channel],'BHZHHZ'), trace = trace(2); end
        if strcmp([trace.channel],'HHZBHZ'), trace = trace(1); end
        if strcmp([trace.channel],'BHEBHNBHZHHEHHNHHZ'), trace = trace(4:6); end
        if strcmp([trace.channel],'HHEHHNHHZBHEBHNBHZ'), trace = trace(1:3); end
        
        % station details from trace
        sta_dts = struct('name',trace(1).station,...
                         'slat',trace(1).latitude,...
                         'slon',trace(1).longitude,...
                         'selev',trace(1).elevation);
        
        % channels' details
        comp = cell(1,length(trace)); for ic = 1:length(trace), comp{ic} = trace(ic).channel(end); end
        chan_dts = struct('name',{{trace.channel}},...
                          'component',{comp},...
                          'sensitivity',[trace.sensitivity],...
                          'dip',[trace.dip],...
                          'azimuth',[trace.azimuth]);
                      
        [gcarc,esaz] = distance(elat,elon,sta_dts.slat,sta_dts.slon);
        [~,seaz] = distance(sta_dts.slat,sta_dts.slon,elat,elon);   

        % DETREND
        for id = 1:length(trace)
            trace(id).data = detrend(trace(id).data);
        end
        
        % SAMPRATE
        if resamprate
%            fprintf('need to apply a low-pass filter to the data too!\n')
            for id = 1:length(trace)
                % anti-aliasing filter - low pass below resamprate/2             
                trace(id).data = filt_quick( trace(id).data,...
                trace(id).sampleRate./trace(id).sampleCount,...
                resamprate/2,1./trace(id).sampleRate);
            end
            samprate = resamprate;
        else
            samprate = round(unique([trace.sampleRate]));
            if length(samprate)>1 
                fprintf(' different samprates, downsamp to min'); 
                samprate = round(min(unique([trace.sampleRate])));
            end
        end
        
      
        % TIME
        tt0 = STARTtime + 1/spd;
        tt1 = ENDtime - 1/spd - 1./samprate/spd;
        tt = [tt0:(1./samprate/spd):tt1]'; 
        
        % SAFETY
        if any([trace.startTime]>tt0)
            fprintf('REQUESTED DATA ONLY STARTS AFTER DESIRED WINDOW START!!\n')
            continue            
        end
        if any([trace.endTime]<tt1)
            fprintf('REQUESTED DATA ENDS BEFORE DESIRED WINDOW END!!\n')
            continue
        end
        
        % data
        nsamps = length(tt);
        dat = zeros(nsamps,length(trace));
        for id = 1:length(trace)
        dat(:,id) = interp1(linspace(trace(id).startTime,...
                                     trace(id).endTime,...
                                     trace(id).sampleCount),...
                                     trace(id).data,    tt);
        end
               
        % fix nans
        nandat = find(isnan(dat));
        if ~isempty(nandat)
            if length(nandat) > 2*length(trace)
                fprintf('lots of nans - look out')
            end
            fprintf('fixing nans')
            dat(nandat) = 0;
        end
            
        
        % phase and distance details
        TT = tauptime('event',[elat,elon],'depth',edep,'station',[sta_dts.slat,sta_dts.slon],'phases',phases);
        for ip = 1:length(TT)
            TT(ip).artime = TT(ip).time/spd + evtime;
        end
        
        % finally put into structure
        data = struct('station',sta_dts,'network',trace(1).network,...
                       'chans',chan_dts,'samprate',samprate,'nsamps',nsamps,...
                       'gcarc',gcarc,'seaz',seaz,'esaz',esaz,...
                       'phases',TT, 'dat',dat,'tt',tt,...
                       'raw',struct('chans',chan_dts,'dat',dat),...
                       'NEZ',false,'rmresp',false,'rmtilt',false,'rmcomp',false);

        nstas4evt = length(datinfo);
        datinfo(nstas4evt+1,1) = struct('sta',stainfo.stas{is},'nwk',stainfo.nwk{is},'chans',{chan_dts.component},'NEZ',false,'rmresp',false,'rmtilt',false,'rmcomp',false,'spectra',false);
        % save
        save(datafile,'data')
        save(datinfofile,'datinfo')
        fprintf('\n')
    end %loop on stas
    kill = []; for ii = 1:length(datinfo), if isempty(datinfo(ii).sta), kill = [kill;ii]; end; end
    datinfo(kill) = [];
    save(datinfofile,'datinfo')
end % loop on evts
cd(wd);

