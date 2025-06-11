%% Script to add new events and stations to existing database 

dbname = 'EARdb';
dbpath = '~/Dropbox/Work/EARdb/';
overwriteall = false;

%% Station parameters
sta_latlims = [-2 17]; % [min_lat max_lat] for stations
sta_lonlims = [32 47]; % [min_lon max_lon] for stations
sta_chans = 'BH*,HH*'; % channel codes to search for

%% Event parameters
mag_lims = [6.0 8.0];
dep_lims = [0 1000]; % set to [0 1000] by default (km)
cpt      = [7.5 40]; % [lat lon] of centre of array;
gc_lims  = [20 135];
startafter = '1990-01-01 00:00:00'; % earliest evtime (yyyy-mm-dd HH:MM:SS)

%% ID for IRIS DMC request
IRIS_ID = 'zeilon';

%% GET TO WORK
wd = pwd;
addpath('matguts');

%% Name directory structure and get all in place
% main database directory
dbpath = regexprep(dbpath,'~',getenv('HOME'));
if ~strcmp(dbpath(end),'/'),dbpath = [dbpath,'/']; end
dbdir = [dbpath];
% info files directory 
infodir = [dbdir,'INFO/'];
% response files directory 
respdir = [dbdir,'INFO/RESP/'];
% data files directory 
datadir = [dbdir,'DATA/'];
% results files directory 
resdir = [dbdir,'RESULTS/'];

cd(dbdir)
%% run project startup file
dbstartup = [dbname,'_startup.m'];
run(dbstartup);
return 

%% Get event information
events_IRIS_new = irisFetch.Events('minimumMagnitude',mag_lims(1),'maximumMagnitude',mag_lims(2),...
                          'minimumDepth',dep_lims(1),'maximumDepth',dep_lims(2),...
                          'startTime',startafter,...
                          'radialcoordinates',[cpt,fliplr(gc_lims)]);
events_request_new = struct('mag_lims',mag_lims,'dep_lims',dep_lims,'startafter',startafter,'centre_point',cpt,'gc_lims',gc_lims);                      

% load existing evinfo
load([infodir,'/events'])

% find new events
[~,ienew] = setdiff(datenum({events_IRIS_new.PreferredTime}'),evinfo.evtimes);
if isempty(ienew), fprintf('0 new events found\n'); return; end
% some info
fprintf('%.0f new events found\n',length(ienew));
ienew_recent = ienew(datenum({events_IRIS_new(ienew).PreferredTime}')>max(evinfo.evtimes));
fprintf('of which \n%.0f are more recent than anything in evinfo\n',length(ienew_recent))

% query whether or not to add new events
ans = input('Add in all (a), recent (r), or NO (anything else) new events? ','s');
if strcmp(ans,'a')
    ieadd = ienew;
elseif strcmp(ans,'r')
    ieadd = ienew_recent;
else
    ieadd = [];
end
ieadd = flipud(ieadd); % make newest to oldest

% add in new events
norids_current = evinfo.norids;
norids_new = norids_current + length(ieadd);
evinfo = struct('orids',[evinfo.orids; norids_current + [1:length(ieadd)]'],...
                'norids',norids_new,...
                'elats',[evinfo.elats;[events_IRIS_new(ieadd).PreferredLatitude]'],...
                'elons',[evinfo.elons;[events_IRIS_new(ieadd).PreferredLongitude]'],...
                'edeps',[evinfo.edeps;[events_IRIS_new(ieadd).PreferredDepth]'],...
                'evmags',[evinfo.evmags;[events_IRIS_new(ieadd).PreferredMagnitudeValue]'],...
                'evtimes',[evinfo.evtimes;datenum({events_IRIS_new(ieadd).PreferredTime}')],...
                'evtimes_IRISstr',{cat(1,evinfo.evtimes_IRISstr,{events_IRIS_new(ieadd).PreferredTime}')},...                 
                'datestamp',cat(1,evinfo.datestamp,datestr(datenum({events_IRIS_new(ieadd).PreferredTime}'),'yyyymmddHHMM')));          

% decide to overwrite
ans = input('Confirm (y) you want to overwrite and add these events? ','s');
if strcmp(ans,'y')
    nupdates = length(events_request);
    events_IRIS = events_IRIS_new;       
    events_request(nupdates+1) = events_request_new;
    % just in case, save old one
    movefile([infodir,'/events.mat'],[infodir,'/events_old.mat']);
    % save new one
    save([infodir,'/events'],'evinfo','events_IRIS','events_request');
end

return
% %% ALTERNATIVELY make from data
% evinfo = mk_evinfo_from_data(datadir);
% save([infodir,'/events'],'evinfo','events_IRIS','events_request');
% return
%% Get station + channel information
stations_IRIS_new = irisFetch.Stations('CHANNEL','*','*','*',sta_chans,...
                            'MinimumLatitude',sta_latlims(1),...
                            'MaximumLatitude',sta_latlims(2),...
                            'MinimumLongitude',sta_lonlims(1),...
                            'MaximumLongitude',sta_lonlims(2));
                        
% parse ongoing stations
for is = 1:length(stations_IRIS_new)
    if isempty(stations_IRIS_new(is).EndDate)
        stations_IRIS_new(is).EndDate = '2599-12-31 23:59:59.000';
    end
end


stations_request_new = struct('lat_lims',sta_latlims,'lon_lims',sta_lonlims,'chans',sta_chans);

stainfo_new = struct('stas',{{stations_IRIS_new.StationCode}'},...
                 'nwk',{{stations_IRIS_new.NetworkCode}'},...
                 'slats',[stations_IRIS_new.Latitude]',...
                 'slons',[stations_IRIS_new.Longitude]',...
                 'selevs',[stations_IRIS_new.Elevation]',...
                 'ondate', datenum({stations_IRIS_new.StartDate}'),...
                 'offdate',datenum({stations_IRIS_new.EndDate}'),...
                 'ondate_str',{{stations_IRIS_new.StartDate}'},...
                 'offdate_str',{{stations_IRIS_new.EndDate}'},...
                 'nstas',length(stations_IRIS_new));  
             

% parse channels             
chans = cell(stainfo_new.nstas,3);
chandips = nan(stainfo_new.nstas,3);
chanazs = nan(stainfo_new.nstas,3);
nchans = zeros(stainfo_new.nstas,1);
for is = 1:stainfo_new.nstas
    nchan = length(stations_IRIS_new(is).Channels);
    tempchans = cell(1,nchan);
    tempdips = nan(1,nchan);
    tempazs = nan(1,nchan);
    for ic = 1:nchan
        tempchans(ic) = {stations_IRIS_new(is).Channels(ic).ChannelCode};
        tempdips(ic) = stations_IRIS_new(is).Channels(ic).Dip;
        tempazs(ic) = stations_IRIS_new(is).Channels(ic).Azimuth;
    end
    [stachans,indch] = unique(tempchans);
    nchans(is) = length(stachans);
    chans(is,1:nchans(is)) = stachans;
    chandips(is,1:nchans(is)) = tempdips(indch);
    chanazs(is,1:nchans(is)) = tempazs(indch);
end
chandips(cellfun('isempty',chans)) = nan;
chanazs(cellfun('isempty',chans)) = nan;
stainfo_new.nchans = nchans;
stainfo_new.chans = chans;
stainfo_new.chandips = chandips;
stainfo_new.chanazs = chanazs;

[stainfo_new] = stainfo_unique(stainfo_new);

% Make plots
figure(44),clf, hold on
m_proj('lambert','lon',sta_lonlims,'lat',sta_latlims); 
[CS,CH]=m_etopo2('contourf',[-5000:500:0 250:250:3000],'edgecolor','none');
 m_grid('linestyle','none','tickdir','out','linewidth',3);
m_coast('color',[0 .2 0],'linewidth',2.5);
colormap([ m_colmap('blues',80); m_colmap('gland',48)]);
brighten(.5);
ax=m_contfbar(1,[.5 .8],CS,CH);
m_plot(stainfo_new.slons,stainfo_new.slats,'or','markerfacecolor','r')
m_plot(stainfo.slons,stainfo.slats,'ob','markerfacecolor','b')




% concat new with existing, find all unique
load([infodir,'/stations']);
stafns = fieldnames(stainfo);
for iff = 1:length(stafns)
    if isscalar(stainfo.(stafns{iff})), continue,
    elseif iscell(stainfo.(stafns{iff}))
        stainfo.(stafns{iff}) = cat(1,stainfo.(stafns{iff}),stainfo_new.(stafns{iff}));
    else
        stainfo.(stafns{iff}) = cat(1,stainfo.(stafns{iff}),stainfo_new.(stafns{iff}));
    end
end
[stainfo] = stainfo_unique(stainfo);



                        
% decide to overwrite
ans = input('Confirm (y) you want to overwrite and add these stations? ','s');
if strcmp(ans,'y')
    nupdates = length(stations_request);
    stations_IRIS = stations_IRIS_new;       
    stations_request(nupdates+1) = stations_request_new;
    % just in case, save old one
    movefile([infodir,'/stations.mat'],[infodir,'/stations_old.mat']);
    % save new one
    save([infodir,'/stations'],'stainfo','stations_IRIS','stations_request');
end

return
%% Append new response files, update all reasponse files
load([infodir,'/stations'],'stainfo');

% grab responses from IRIS
stations_IRIS_r = irisFetch.Stations('RESPONSE','*','*','*',sta_chans,...
                            'MinimumLatitude',sta_latlims(1),...
                            'MaximumLatitude',sta_latlims(2),...
                            'MinimumLongitude',sta_lonlims(1),...
                            'MaximumLongitude',sta_lonlims(2));


[stainfo] = stainfo_extractfromIRIS_resp(stations_IRIS_r,stainfo)

% decide to overwrite
ans = input('Confirm (y) you want to overwrite and add these responses? ','s');
if strcmp(ans,'y')
    % just in case, save old one
    copyfile([infodir,'/stations.mat'],[infodir,'/stations_noresp.mat']);
    % save new one
    save([infodir,'/stations'],'stainfo','-append');
end

return 

%% Response SAC_PZ files
% Build and send BREQFAST request file for dataless seed                     
addpath('~/Dropbox/MATLAB/seis_tools/breqfasting/');
load([infodir,'/stations'])
breq_fast_request([dbname,'_dataless'],IRIS_ID,stainfo.stas,{'*'},stainfo.nwk,'',stainfo.ondate_str,stainfo.offdate_str,'dataless_SEED',dbname);
return 
%% ======= WAIT FOR NOTIFICATION THAT DATALESS SEED IS ON SERVER ======= %%
% Download and process dataless seed  
breq_fast_dataless_PZprocess([dbname,'_dataless'],IRIS_ID,respdir,{'BH*','HH*'},1 )     
return
                        