%% Script to get some details about the database we just built

% project details
dbname = 'EARdb';
dbdir = '/Users/zeilon/Work/EastAfrica/EARdb/'; % include final slash

dbname = 'FRES_PILOT';
dbdir = '~/Dropbox/Work/FRES_PILOT/'; % include final slash


%% Preliminaries
wd = pwd;
addpath('matguts')
cd(dbdir);
run([dbdir,dbname,'_startup.m']);
addpath('~/Dropbox/MATLAB/lib/m_map/');
run('~/Dropbox/MATLAB/lib/seizmo/startup_seizmo.m');
% event details
load([infodir,'/events'],'evinfo'); 
% station details
load([infodir,'/stations'],'stainfo');

%% Histogram of station recordings
t0 = min(stainfo.ondate);
t1 = ceil(now);
tt = t0:t1;
Nonstas = zeros(length(tt),1);
for it = 1:length(tt)
    Nonstas(it) = sum(stainfo.ondate<tt(it) & stainfo.offdate>tt(it));
end
figure(22)
bar(tt,Nonstas)
datetick('x')
    
%% Polar plot of earthquakes
figure(23), clf
m_proj('stereographic','lat',mean(stainfo.slats),'long',mean(stainfo.slons),'radius',115);
% m_elev('contour',[-3500:1000:-500],'edgecolor','b');
m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_range_ring(mean(stainfo.slons),mean(stainfo.slats),[1100:1100:22000],'color',[.3 .8 .5],'linewi',1.5);
m_line(evinfo.elons,evinfo.elats,'marker','o','color',[0 0 .5],...
          'linest','none','markerfacecolor','none','clip','point');
m_grid('xtick',[],'tickdir','out','ytick',[],'linest','-','linewidth',2);

%% Area plot of Stations
figure(24),clf
latlims = [min(stainfo.slats)-1, max(stainfo.slats)+1];
lonlims = [min(stainfo.slons)-1, max(stainfo.slons)+1];
m_proj('Gall-Peters','lat',latlims,'long',lonlims);
m_coast('patch',[.7 .7 .7],'edgecolor','none');
[~,he] = m_elev; set(he,'linestyle','--','linewidth',1.1,'color',0.2*[1 1 1])
m_gshhs_i('color','k')

m_line(stainfo.slons,stainfo.slats,'marker','^','color','r',...
          'linest','none','markerfacecolor','r','markeredgecolor','k','clip','point','linewidth',0.4);
m_grid('linest','-','linewidth',2);
m_etopo2


