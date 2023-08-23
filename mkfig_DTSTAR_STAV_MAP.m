% script to plot mean differential attenuation at each station. The mean
% value at each station is solved for in a least squares sense by
% considering all individual measurements at each station for each event
% and accounting for event terms in the function lsq_sta_evt - you may want
% to just take a simple arithmetic average of the station's individual
% measurements.
% 
% Z. Eilon 2016


clear all
close all
addpath('matguts')
% project details
dbname = 'EARdb';
dbdir = '/Users/zeilon/Work/EastAfrica/EARdb/'; % include final slash



%% parameters
ifsave = false;

method = 'comb';% 'comb' or 'specR'

ifOBSinv = false;
ifOBSonly = true;

plotsize = 800;
scale = 100; % length of lines
phases = {'P','S'};
components = {'Z','T'}; %'Z', 'R', or 'T'

ALP = []; % [] for the default (also 0), or the value of alpha

keyloc = [-131.1,49.1];
keysiz = [1.3,1.8];

%% directories 
% ANTELOPE DB DETAILS
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% DATA DIRECTORY (top level)
datadir = '/Volumes/DATA/CASCADIA/DATA/'; % needs final slash
% RESULTS DIRECTORY
resdir = '~/Documents/MATLAB/CASC_atten/results/'; % needs final slash
% FIGURES DIRECTORY
figdir = 'figs';

cmap = parula;
% tstlim = 2.*[-1.5 0.5];

close all
%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

obsloadstr = ''; 
if ifOBSinv, obsloadstr = 'OBS_'; obsstr = 'OBSinv_'; end
 
if isempty(ALP), alpstr = ''; else alpstr = sprintf('_alp%03.0f',ALP*100); end

%% GET EVENTS+STATIONS DATA
[ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam );
[ nstas_all,stas_all,slats_all,slons_all,selevs_all,~,~,~,statypes_all ] = db_stadata( dbdir,dbnam );



%% -----------  LOOP THROUGH PHASES/COMPS  ------------
for ip = 1:length(phases)
phase = phases{ip};
component = components{ip};
tstlim = ip*[-0.8 0.8];

%% LOAD RESULTS
if strcmp(method,'specR')
    load([resdir,'all_dtstar_',obsloadstr,phase,'_',component,alpstr,'.mat']); % alpstr has to be '' for this - just a check.
elseif strcmp(method,'comb')
    load([resdir,'all_dtstar',method,'_',obsloadstr,phase,'_',component,alpstr,'.mat']);
    all_dtstar = all_dtstar_comb;
else
    error('Need to specify method')
end

%% limit to only obs if needed
if ifOBSonly
    all_dtstar(~strcmp(statypes_all,'OBS'),:) = nan;
    obsstr = 'OBSonly_';
else 
    obsstr = '';
end
    

%% parse the stations and their data - collate repeats and remove allnan stas
[ stas,all_dtstar,slats,slons,selevs ] = results_PARSE_STATIONS( stas_all,all_dtstar,slats_all,slons_all,selevs_all );
nstas = length(stas);
sages = jdf_crust_age( slats,slons );

%% compute station averages
% nnan = ~isnan(nanmean(all_dtstar,2));
Nobs = sum(~isnan(all_dtstar),2);
% compute station averages
[ sta_terms,evt_terms ] = lsq_sta_evt( all_dtstar,0.01);
% figure(4), hist(evt_terms)

%% ------------------ MAP WITH DT FOR THIS EVENT ------------------

figure(31), clf, hold on
mkfig_CascMAP
lonlims = [-132.1 -120];
set(gcf,'position',[200 300 plotsize/plot_size_ratio(lonlims,latlims) plotsize])
set(gca,'xlim',lonlims)

% nnan = ~isnan(nanmean(all_dtstar,2));
Nobs = sum(~isnan(all_dtstar),2);
% compute station averages
[ sta_terms,evt_terms ] = lsq_sta_evt( all_dtstar,0.01);

%% plot the station-averaged differential tstar
hp = scatter(slons,slats,scale*sqrt(Nobs),sta_terms,'filled','MarkerEdgeColor','k');

colormap(cmap), caxis(tstlim)

%% colour bar
tkvl = unique(round_level([tstlim(1):0.5:tstlim(2)],0.5));
cbar_custom(gca, 'location',[-131.4 -131. 39.7 43],'tickside','right',...
    'lims',tstlim,'tickvals',tkvl,'cmap',cmap,...
    'FontSize',12,'FontWeight','bold',...
	'title',sprintf('$\\Delta t^*_%s$ \\,(s)',phase),'interpreter','latex');

%% key
kle = keyloc(1) - 0.5*keysiz(1); kri = keyloc(1) + 0.5*keysiz(1);
kbo = keyloc(2) - 0.5*keysiz(2); kto = keyloc(2) + 0.5*keysiz(2);
kw = keysiz(1); kh = keysiz(2);

hold on
patch([kle kle kri kri kle],[kbo kto kto kbo kbo],...
      'w','LineWidth',2)
text(keyloc(1),kto-0.1*kh,'Nevts',...
    'fontsize',15,'fontweight','bold','verticalalignment','top','horizontalalignment','center');

scatter(kle + 0.3*kw,kbo + 0.61*kh,scale*sqrt(1),  'k','filled','MarkerEdgeColor','k');
scatter(kle + 0.3*kw,kbo + 0.435*kh,scale*sqrt(10), 'k','filled','MarkerEdgeColor','k');
scatter(kle + 0.3*kw,kbo + 0.18*kh,scale*sqrt(50),'k','filled','MarkerEdgeColor','k');

text(kle + 0.62*kw,kbo + 0.61*kh,'1', 'fontsize',12,'fontweight','bold','verticalalignment','middle');
text(kle + 0.62*kw,kbo + 0.435*kh,'10','fontsize',12,'fontweight','bold','verticalalignment','middle');
text(kle + 0.62*kw,kbo + 0.18*kh,'50','fontsize',12,'fontweight','bold','verticalalignment','middle');

%% title
title(sprintf('%s $\\Delta t^*$ for $%s$-waves (%s component) %s',strtok(obsstr,'_'),phase,component,method),...
      'FontSize',18,'FontWeight','bold','Interpreter','latex')
set(gca,'FontSize',14,'LineWidth',2.5,'box','on')

% save
if ifsave
save2pdf(31,sprintf('MAP_dtstar_staav_%s_%s%s_%s%s',method,obsstr,phase,component,alpstr),figdir);

results = struct('stas',{stas},'dtstar',sta_terms,'slats',slats,'slons',slons,'selevs',selevs,'Nobs',Nobs);
resfile = sprintf('stav_dtstar%s_%s%s_%s%s',method,obsstr,phase,component,alpstr);
save([resdir,resfile],'results')
end
return
end % loop on 2 phases