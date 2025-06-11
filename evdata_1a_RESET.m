% script to go through data and reset to raw traces and channels=
clear

overwrite = true;

% project details
dbname = 'EARdb';
dbdir = '/Users/zeilon/Dropbox/Work/EARdb/'; % include final slash

%% Preliminaries
wd = pwd;
addpath('matguts')
cd(dbdir);
run([dbdir,dbname,'_startup.m']);
javaaddpath('IRIS-WS-2.0.15.jar')
spd = 24*60*60; % seconds per day - to convert time in seconds to serial time

% event details
load([infodir,'/events'],'evinfo'); 
% station details
load([infodir,'/stations'],'stainfo');

for ie = 3651%:evinfo.norids % 1:norids % loop on orids
    orid = evinfo.orids(ie);
    fprintf('\n Orid %.0f %s \n\n',orid,evinfo.evtimes_IRISstr{ie}(1:end-4))
    evdir = [num2str(orid,'%03d'),'_',char(evinfo.datestamp(ie,:)),'/'];
    datinfofile = [datadir,evdir,'_datinfo'];
   
    % check files exist
    if ~exist([datinfofile,'.mat'],'file'), fprintf('No station mat files for this event\n');continue, end
    load(datinfofile)
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end

    for is = 1:length(datinfo) % loop on stas
        sta = datinfo(is).sta; % sta name        
        load([datadir,evdir,sta,'.mat']); % load sta data for this evt
        fprintf('Station %.0f %-5s...',is,sta)
        data.dat    = data.raw.dat;
        
        % fix nans
        nandat = find(isnan(data.dat));
        if ~isempty(nandat)
            if length(nandat) > 2*size(data.dat,2)
                fprintf('lots of nans - look out ')
            end
            fprintf('fixing nans')
            data.dat(nandat) = 0;
        end
            
        data.chans  = data.raw.chans;
        data.NEZ    = false;
        data.rmresp = false;
        data.rmtilt = false;
        data.rmcomp = false;
        
        [data.phases.xtime]   = deal([]);
        [data.phases.xartime] = deal([]);
        [data.phases.xacor]   = deal([]);
        
        datinfo(is).chans  = data.raw.chans.component;
        datinfo(is).NEZ    = false;
        datinfo(is).rmresp = false;
        datinfo(is).rmtilt = false;
        datinfo(is).rmcomp = false;
        
        %% log reset and save
        fprintf(' reset\n')
        save([datadir,evdir,sta],'data')
        save(datinfofile,'datinfo')
        
        copyfile([datinfofile,'.mat'],[datinfofile,'_P.mat'])
        copyfile([datinfofile,'.mat'],[datinfofile,'_S.mat'])
            
    end % loop on stas

	%% sum up
    fprintf(' STA  CHAN  NEZ  resp  spectra  tilt  comp\n')
    for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f      %1.0f       %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).spectra,datinfo(is).rmtilt,datinfo(is).rmcomp); end
    
end
cd(wd);
