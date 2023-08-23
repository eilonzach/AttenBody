%% Script to establish a database of stations and events for body wave study

dbname = 'FRES_PILOT';
dbpath = '~/Dropbox/Work/'; % no need for name of db here, this is where it will sit

%% Station parameters
sta_latlims = [37 41]; % [min_lat max_lat] for stations
sta_lonlims = [-122.5 -119]; % [min_lon max_lon] for stations
sta_chans = 'BH*,HH*'; % channel codes to search for

%% Event parameters
mag_lims = [6.0 8.0];
dep_lims = [0 1000]; % set to [0 1000] by default (km)
cpt      = [7.5 40]; % [lat lon] of centre of array;
gc_lims  = [20 135];
startafter = '2000-01-01 00:00:00'; % earliest evtime (yyyy-mm-dd HH:MM:SS)

%% ID for IRIS DMC request
IRIS_ID = 'zeilon';


%% GET TO WORK
wd = pwd;
addpath('matguts');

%% Make directory structure
% main database directory
dbpath = regexprep(dbpath,'~',getenv('HOME'));
if ~strcmp(dbpath(end),'/'),dbpath = [dbpath,'/']; end
dbdir = [dbpath,dbname,'/'];
if exist(dbdir,'dir')~=7, mkdir(dbdir); end
% info files directory 
infodir = [dbdir,'INFO/'];
if exist(infodir,'dir')~=7, mkdir(infodir); end
% response files directory 
respdir = [dbdir,'INFO/RESP/'];
if exist(respdir,'dir')~=7, mkdir(respdir); end
% data files directory 
datadir = [dbdir,'DATA/'];
if exist(datadir,'dir')~=7, mkdir(datadir); end
% results files directory 
resdir = [dbdir,'RESULTS/'];
if exist(resdir,'dir')~=7, mkdir(resdir); end

cd(dbdir)

%% make project startup file
dbstartup = [dbname,'_startup.m'];
fid = fopen(dbstartup,'w');
fprintf(fid,'%%%% Startup file for project %s\n',dbname);
fprintf(fid,'global dbdir infodir respdir datadir resdir;\n');
fprintf(fid,'dbdir = ''%s'';\n',dbdir);
fprintf(fid,'infodir = ''%s'';\n',infodir);
fprintf(fid,'respdir = ''%s'';\n',respdir);
fprintf(fid,'datadir = ''%s'';\n',datadir);
fprintf(fid,'resdir = ''%s'';\n',resdir);
fclose(fid);
run(dbstartup);
return 

addpath('~/Dropbox/MATLAB/seis_tools/')

%% Get event information
fprintf('%s\nDOWNLOADING event information using irisFetch\n',repmat('*',1,60))
events_IRIS = irisFetch.Events('minimumMagnitude',mag_lims(1),'maximumMagnitude',mag_lims(2),...
                          'minimumDepth',dep_lims(1),'maximumDepth',dep_lims(2),...
                          'startTime',startafter,...
                          'radialcoordinates',[cpt,fliplr(gc_lims)]);

events_request = struct('mag_lims',mag_lims,'dep_lims',dep_lims,'startafter',startafter,'centre_point',cpt,'gc_lims',gc_lims);

% invert order so earlier events first
events_IRIS = fliplr(flipud(events_IRIS));

% put in evinfo structure
evinfo = struct('orids',[1:length(events_IRIS)]',...
                'norids',length(events_IRIS),...
                'elats',[events_IRIS.PreferredLatitude]',...
                'elons',[events_IRIS.PreferredLongitude]',...
                'edeps',[events_IRIS.PreferredDepth]',...
                'evmags',[events_IRIS.PreferredMagnitudeValue]',...
                'evtimes',datenum({events_IRIS.PreferredTime}'),...
                'evtimes_IRISstr',{{events_IRIS.PreferredTime}'},...                 
                'datestamp',datestr(datenum({events_IRIS.PreferredTime}'),'yyyymmddHHMM'));                 
%evinfo_mag_dist_cull
save([infodir,'/events'],'evinfo','events_IRIS','events_request');

return
%% ALTERNATIVELY make from data
evinfo = mk_evinfo_from_data(datadir);
save([infodir,'/events'],'evinfo','events_IRIS','events_request');
return
%% Get station + channel information
fprintf('%s\nDOWNLOADING station information using irisFetch\n',repmat('*',1,60))

addpath('~/Dropbox/MATLAB/seis_tools/')

stations_request = struct('lat_lims',sta_latlims,'lon_lims',sta_lonlims,'chans',sta_chans);

% stations_IRIS = irisFetch.Stations('CHANNEL','*','*','*',sta_chans,...
%                             'MinimumLatitude',sta_latlims(1),...
%                             'MaximumLatitude',sta_latlims(2),...
%                             'MinimumLongitude',sta_lonlims(1),...
%                             'MaximumLongitude',sta_lonlims(2));
% [stainfo] = stainfo_extractfromIRIS(stations_IRIS,stainfo)


stations_IRIS_r = irisFetch.Stations('RESPONSE','*','*','*',sta_chans,...
                            'MinimumLatitude',sta_latlims(1),...
                            'MaximumLatitude',sta_latlims(2),...
                            'MinimumLongitude',sta_lonlims(1),...
                            'MaximumLongitude',sta_lonlims(2));


[stainfo] = stainfo_extractfromIRIS_resp(stations_IRIS_r,stainfo)

% % find and fill in empty enddates
% for is = 1:length(stations_IRIS)
%     if isempty(stations_IRIS(is).EndDate)
%         stations_IRIS(is).EndDate = '2599-12-31 23:59:59.999';
%     end
% end
% 
% 
% stainfo = struct('stas',{{stations_IRIS.StationCode}'},...
%                  'nwk',{{stations_IRIS.NetworkCode}'},...
%                  'slats',[stations_IRIS.Latitude]',...
%                  'slons',[stations_IRIS.Longitude]',...
%                  'selevs',[stations_IRIS.Elevation]',...
%                  'ondate',datenum({stations_IRIS.StartDate}'),...
%                  'offdate',datenum({stations_IRIS.EndDate}'),...
%                  'ondate_str',{{stations_IRIS.StartDate}'},...
%                  'offdate_str',{{stations_IRIS.EndDate}'},...
%                  'nstas',length(stations_IRIS));  
%              
% % parse channels             
% chans = cell(stainfo.nstas,3);
% chandips = nan(stainfo.nstas,3);
% chanazs = nan(stainfo.nstas,3);
% nchans = zeros(stainfo.nstas,1);
% for is = 1:stainfo.nstas
%     nchan = length(stations_IRIS(is).Channels);
%     tempchans = cell(1,nchan);
%     tempdips = nan(1,nchan);
%     tempazs = nan(1,nchan);
%     for ic = 1:nchan
%         tempchans(ic) = {stations_IRIS(is).Channels(ic).ChannelCode};
%         tempdips(ic) = stations_IRIS(is).Channels(ic).Dip;
%         tempazs(ic) = stations_IRIS(is).Channels(ic).Azimuth;
%     end
%     [stachans,indch] = unique(tempchans);
%     nchans(is) = length(stachans);
%     chans(is,1:nchans(is)) = stachans;
%     chandips(is,1:nchans(is)) = tempdips(indch);
%     chanazs(is,1:nchans(is)) = tempazs(indch);
% end
% chandips(cellfun('isempty',chans)) = nan;
% chanazs(cellfun('isempty',chans)) = nan;
% stainfo.nchans = nchans;
% stainfo.chans = chans;
% stainfo.chandips = chandips;
% stainfo.chanazs = chanazs;
% 
% [stainfo] = stainfo_unique(stainfo);

                        
save([infodir,'/stations'],'stainfo','stations_IRIS','stations_request');

%% distances
fprintf('Computing all inter-station distances\n');
dd = zeros(stainfo.nstas);
for is1 = 1:stainfo.nstas
for is2 = is1+1:stainfo.nstas
    dd(is1,is2) = distance(stainfo.slats(is1),stainfo.slons(is1),stainfo.slats(is2),stainfo.slons(is2));
end
end
dd = dd+dd';
save([infodir,'/stations'],'stainfo','dd','stations_IRIS','stations_request');


return

%% Response SAC_PZ files
% Build and send BREQFAST request file for dataless seed                     
addpath('~/Dropbox/MATLAB/seis_tools/breqfasting/');
breq_fast_request([dbname,'_dataless'],IRIS_ID,{stations_IRIS.StationCode}',{'*'},{stations_IRIS.NetworkCode}','',{stations_IRIS.StartDate}',{stations_IRIS.EndDate}','dataless_SEED',dbname);
return 
%% ======= WAIT FOR NOTIFICATION THAT DATALESS SEED IS ON SERVER ======= %%
% Download and process dataless seed  
breq_fast_dataless_PZprocess([dbname,'_dataless'],IRIS_ID,respdir,{'BH*','HH*'},1 )     
return
                        