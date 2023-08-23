clear all
close all
cd ~/Documents/MATLAB/CASC_atten/
addpath('matguts')

%% PARMS
ifsave = false;
specmethod = 'comb';% 'comb' or 'specR'
indiv_or_stav = 'stav'; % if stav, load the results structures made in the mkfig_DT/TSTAR_MAP scripts
ifOBSonly = true; % to get only OBS results


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

% station details
[ nstas,stas,slats,slons,selevs,ondate,offdate,staname,statype ] = db_stadata( dbdir,dbnam );

%% =================================================================== %%
%% ==========================  INDIVIDUAL  =========================== %%
%% =================================================================== %%
if strcmp(indiv_or_stav,'indiv')
% LOAD INDIVIDUAL RESULTS
load([resdir,'all_dT_P_Z.mat']); dTp = all_dT;
load([resdir,'all_dT_S_T.mat']); dTs = all_dT;
if strcmp(specmethod,'specR')
    load([resdir,'all_dtstar_P_Z.mat']); dtstp = all_dtstar;
    load([resdir,'all_dtstar_S_T.mat']); dtsts = all_dtstar;
elseif strcmp(specmethod,'comb')
    load([resdir,'all_dtstarcomb_P_Z.mat']); dtstp = all_dtstar_comb;
    load([resdir,'all_dtstarcomb_S_T.mat']); dtsts = all_dtstar_comb;
else
    error('Need to specify method of obtaining dtstar')
end

norids = size(all_dT,2);
isob = (strcmp(statype,'OBS' )*ones(1,norids))==1;
isla = (strcmp(statype,'LAND')*ones(1,norids))==1;


%% =====================  INDIVIDUAL MEASUREMENTS  ===================== %%
ydT  = dTp~=0   & dTs~=0   & ~isnan(dTp)   & ~isnan(dTs)   & abs(dTp)<3   & abs(dTs)<5;
ydts = dtstp~=0 & dtsts~=0 & ~isnan(dtstp) & ~isnan(dtsts) & abs(dtstp)<3 & abs(dtsts)<5;
yP   = dTp~=0 & ~isnan(dTp) & dtstp~=0 & ~isnan(dtstp) & abs(dTp)<3 & abs(dtstp)<3;
yS   = dTs~=0 & ~isnan(dTs) & dtsts~=0 & ~isnan(dtsts) & abs(dTs)<5 & abs(dtsts)<5;

isob = (strcmp(statype,'OBS' )*ones(1,norids))==1;
isla = (strcmp(statype,'LAND')*ones(1,norids))==1;

fprintf('For S: %.0f individual dT measurements\n',sum(sum(dTs~=0 & ~isnan(dTs) &abs(dTs)<5 )))
fprintf('For P: %.0f individual dT measurements\n',sum(sum(dTp~=0 & ~isnan(dTp) &abs(dTp)<3 )))
fprintf('... giving %.0f combined P-S dT pairs\n',sum(sum(ydT)));
fprintf('For S: %.0f individual dt* measurements\n',sum(sum(dtsts~=0 & ~isnan(dtsts) &abs(dtsts)<5 )))
fprintf('For P: %.0f individual dt* measurements\n',sum(sum(dtstp~=0 & ~isnan(dtstp) &abs(dtstp)<3 )))
fprintf('... giving %.0f combined P-S dt* pairs\n',sum(sum(ydts)));

% diff-tt line fit
tts = dTs(ydT) - mean( dTs(ydT) );
ttp = dTp(ydT) - mean( dTp(ydT) );
[ m_tt] = fit_LSqEr( ttp,tts,1,2,1:0.01:5,[],2,1);

% diff-tstar line fit
tss = dtsts(ydts) - mean( dtsts(ydts) );
tsp = dtstp(ydts) - mean( dtstp(ydts) );
[ m_ts,b,mstd,bstd ] = fit_LSqEr( tsp,tss,1,1,1:0.01:5,[],4,1);



%% PLOT DIFFERENTIAL TRAVEL TIME
figure(1), clf, set(gcf,'position',[200 200 570 610])
hold on

scatter(ttp,tts,40,'k','o','LineWidth',1.5)
scatter(ttp(isob(ydT)),tts(isob(ydT)),40,'b','o','LineWidth',1.5)

axis([-3.7 3.7 -5 5]); ax = axis;

plot(ax(1:2),m_tt*ax(1:2),'r--','LineWidth',3)
% plot(ax(1:2),2.2*1.81*ax(1:2),'k--','LineWidth',3)
% scatter(ttp,tts,40,'b','o','LineWidth',1.5)

text(-3.2,4.6,sprintf('\\textbf{Slope = %.2f}',m_tt),'FontSize',17,'interpreter','latex')
text(-3,4.15,sprintf('\\textbf{$\\Rightarrow\\, \\delta \\ln{V_S}/\\delta \\ln{V_P} = %.2f$}',m_tt/1.81),'FontSize',17,'interpreter','latex')

set(gca,'fontsize',15,'color','none','LineWidth',2,'box','on')
% title('Ratio of $P$ and $S$ travel times residuals','FontSize',18,'interpreter','latex')
xlabel('$P$ residual time, $\delta T_P$ (s)','FontSize',20,'interpreter','latex')
ylabel('$S$ residual time, $\delta T_S$ (s)','FontSize',20,'interpreter','latex')
if ifsave
save2pdf(1,'dTS_vs_dTP_indiv',figdir); pause(0.1)
end

%% PLOT DIFFERENTIAL T-STAR
figure(2), clf, set(gcf,'position',[300 200 570 610])
hold on

scatter(tsp,tss,40,'k','o','LineWidth',1.5)
scatter(tsp(isob(ydts)),tss(isob(ydts)),40,'b','o','LineWidth',1.5)

axis([-3.7 3.7 -5 5]); ax = axis;

plot(ax(1:2),m_ts*ax(1:2),'r--','LineWidth',3)

text(-3.2,4.6,sprintf('\\textbf{Slope = %.2f}',m_ts),'FontSize',17,'interpreter','latex')
text(-3,4.15,sprintf('\\textbf{$\\Rightarrow\\, Q_P/Q_S = %.2f$}',m_ts/1.81),'FontSize',17,'interpreter','latex')

set(gca,'fontsize',15,'color','none','LineWidth',2,'box','on')
% title('Ratio of $P$ and $S$ travel times residuals','FontSize',18,'interpreter','latex')
xlabel('$P$ differential t-star, $\Delta t^*_P$ (s)','FontSize',20,'interpreter','latex')
ylabel('$S$ differential t-star,  $\Delta t^*_S$ (s)','FontSize',20,'interpreter','latex')
if ifsave
save2pdf(2,'dtstarS_vs_dtstarP_indiv',figdir); pause(0.1)
end

%% PLOT S tt vs. tstar
figure(3), clf, hold on, set(gcf,'position',[300 200 630 570])
scatter(dTs(yS),dtsts(yS),40,'k','o','LineWidth',1.5)
scatter(dTs(isob&yS),dtsts(isob&yS),40,'b','o','LineWidth',1.5)

axis([-4.7 4.7 -4.7 4.7])
set(gca,'fontsize',15,'color','none','LineWidth',2,'box','on')
xlabel('$S$ residual time, $\delta T_S$ (s)','FontSize',20,'interpreter','latex')
ylabel('$S$ differential t-star, $\Delta t^*_S$ (s)','FontSize',20,'interpreter','latex')

if ifsave
save2pdf(3,'dtstar_vs_dT_S_indiv',figdir); pause(0.1)
end

%% PLOT P tt vs. tstar
figure(4), clf, hold on, set(gcf,'position',[300 200 630 570])
scatter(dTp(yP),dtstp(yP),40,'k','o','LineWidth',1.5)
scatter(dTp(isob&yP),dtstp(isob&yP),40,'b','o','LineWidth',1.5)

axis([-2.1 2.1 -2.1 2.1])
set(gca,'fontsize',15,'color','none','LineWidth',2,'box','on')
xlabel('$P$ residual time, $\delta T_P$','FontSize',20,'interpreter','latex')
ylabel('$P$ differential t-star, $\Delta t^*_P$','FontSize',20,'interpreter','latex')

if ifsave
save2pdf(4,'dtstar_vs_dT_P_indiv',figdir); pause(0.1)
end

%% =====================  OBS INDIV. MEASUREMENTS  ===================== %%
dTpo = dTp(strcmp(statype,'OBS'),:);
dTso = dTs(strcmp(statype,'OBS'),:);
dtstpo = dtstp(strcmp(statype,'OBS'),:);
dtstso = dtsts(strcmp(statype,'OBS'),:);

ydTo  = dTpo~=0   & dTso~=0   & ~isnan(dTpo)   & ~isnan(dTso)   & abs(dTpo)<3   & abs(dTso)<5;
ydtso = dtstpo~=0 & dtstso~=0 & ~isnan(dtstpo) & ~isnan(dtstso) & abs(dtstpo)<3 & abs(dtstso)<5;
yPo   = dTpo~=0 & ~isnan(dTpo) & dtstpo~=0 & ~isnan(dtstpo) & abs(dTpo)<3 & abs(dtstpo)<3;
ySo   = dTso~=0 & ~isnan(dTso) & dtstso~=0 & ~isnan(dtstso) & abs(dTso)<5 & abs(dtstso)<5;

% diff-tt line fit
ttso = dTso(ydTo) - mean( dTso(ydTo) );
ttpo = dTpo(ydTo) - mean( dTpo(ydTo) );
[ m_tto] = fit_LSqEr( ttpo,ttso,1,2,1:0.01:5,[],2,1);

% diff-tstar line fit
tsso = dtstso(ydtso) - mean( dtstso(ydtso) );
tspo = dtstpo(ydtso) - mean( dtstpo(ydtso) );
[ m_tso ] = fit_LSqEr( tspo,tsso ,1,1,1:0.01:5,[],2,1);

%% PLOT OBS DIFFERENTIAL TRAVEL TIME
figure(5), clf, set(gcf,'position',[200 200 480 610])
hold on

scatter(ttpo,ttso,40,'b','o','LineWidth',1.5)

axis([-3.7 3.7 -5 5]); ax = axis;

plot(ax(1:2),m_tto*ax(1:2),'r--','LineWidth',3)
% plot(ax(1:2),2.2*1.81*ax(1:2),'k--','LineWidth',3)
% scatter(ttp,tts,40,'b','o','LineWidth',1.5)

text(-3.2,4.6,sprintf('\\textbf{Slope = %.2f}',m_tto),'FontSize',17,'interpreter','latex')
text(-3,4.15,sprintf('\\textbf{$\\Rightarrow\\, \\delta \\ln{V_S}/\\delta \\ln{V_P} = %.2f$}',m_tto/1.81),'FontSize',17,'interpreter','latex')

set(gca,'fontsize',12,'color','none')
% title('Ratio of $P$ and $S$ travel times residuals','FontSize',18,'interpreter','latex')
xlabel('OBS $P$ residual time (s)','FontSize',16,'interpreter','latex')
ylabel('OBS $S$ residual time (s)','FontSize',16,'interpreter','latex')
if ifsave
save2pdf(5,'dTS_vs_dTP_indiv_OBS',figdir); pause(0.1)
end
%% PLOT OBS DIFFERENTIAL T-STAR
figure(6), clf, set(gcf,'position',[300 200 540 610])
hold on

scatter(tspo,tsso,40,'b','o','LineWidth',1.5)

axis([-3.7 3.7 -5 5]); ax = axis;

plot(ax(1:2),m_tso*ax(1:2),'r--','LineWidth',3)

text(-3.2,4.6,sprintf('\\textbf{Slope = %.2f}',m_tso),'FontSize',17,'interpreter','latex')
text(-3,4.15,sprintf('\\textbf{$\\Rightarrow\\, Q_P/Q_S = %.2f$}',m_tso/1.81),'FontSize',17,'interpreter','latex')

set(gca,'fontsize',12,'color','none','LineWidth',2,'box','on')
% title('Ratio of $P$ and $S$ travel times residuals','FontSize',18,'interpreter','latex')
xlabel('OBS $P$-wave $\Delta t^*$ (s)','FontSize',20,'interpreter','latex')
ylabel('OBS $S$-wave $\Delta t^*$ (s)','FontSize',20,'interpreter','latex')

if ifsave
save2pdf(6,'dtstarS_vs_dtstarP_indiv_OBS',figdir); pause(0.1)
end

%% PLOT OBS S tt vs. tstar
figure(7)
scatter(dTso(ySo),dtstso(ySo),40,'b','o')

set(gca,'fontsize',12,'color','none','LineWidth',2,'box','on')
% title('Ratio of differential travel time to t-star for S','FontSize',18,'interpreter','latex')
xlabel('OBS $\delta T_S$','FontSize',15,'interpreter','latex')
ylabel('OBS $\Delta t^*_S$','FontSize',15,'interpreter','latex')

%% PLOT OBS P tt vs. tstar
figure(8)
scatter(dTpo(yPo),dtstpo(yPo),40,'b','o')

set(gca,'fontsize',12,'color','none','LineWidth',2,'box','on')
% title('Ratio of differential travel time to t-star for S','FontSize',18,'interpreter','latex')
xlabel('OBS $\delta T_P$','FontSize',15,'interpreter','latex')
ylabel('OBS $\Delta t^*_P$','FontSize',15,'interpreter','latex')

end % on if indiv

%% =================================================================== %%
%% ========================  STATION AVERAGES  ======================= %%
%% =================================================================== %%
if strcmp(indiv_or_stav,'stav')

% LOAD STAV RESULTS
if ~ifOBSonly
    obsstr = '';
elseif ifOBSonly
    obsstr = '_OBSonly';
end

load([resdir,'stav_dTxcorr',obsstr,'_P_Z.mat']); stavR_dTp = results;
load([resdir,'stav_dTxcorr',obsstr,'_S_T.mat']); stavR_dTs = results;
load([resdir,'stav_dtstar',specmethod,obsstr,'_P_Z.mat']); stavR_dtstarp = results;
load([resdir,'stav_dtstar',specmethod,obsstr,'_S_T.mat']); stavR_dtstars = results;

% parse into comparable vectors, same length
stav_stas = unique([stavR_dtstarp.stas;stavR_dtstars.stas;stavR_dTp.stas;stavR_dTs.stas])';
nstas = length(stav_stas);
for is = 1:nstas
    stav_statype(is) = statype(strcmp(stas,strtok(stav_stas(is),'_')));
end
stav_isob = strcmp(stav_statype,'OBS')';
stav_isla = strcmp(stav_statype,'LAND')';


stav_dTp = nan(nstas,1); Nobs_dTp = nan(nstas,1);
stav_dTs = nan(nstas,1); Nobs_dTs = nan(nstas,1);
stav_dtstarp = nan(nstas,1); Nobs_dtstarp = nan(nstas,1);
stav_dtstars = nan(nstas,1); Nobs_dtstars = nan(nstas,1);
[~,ia,ib] = intersect(stav_stas,stavR_dTp.stas); stav_dTp(ia) = stavR_dTp.dT(ib); Nobs_dTp(ia) = stavR_dTp.Nobs(ib); 
[~,ia,ib] = intersect(stav_stas,stavR_dTs.stas); stav_dTs(ia) = stavR_dTs.dT(ib); Nobs_dTs(ia) = stavR_dTs.Nobs(ib); 
[~,ia,ib] = intersect(stav_stas,stavR_dtstarp.stas); stav_dtstarp(ia) = stavR_dtstarp.dtstar(ib); Nobs_dtstarp(ia) = stavR_dtstarp.Nobs(ib);
[~,ia,ib] = intersect(stav_stas,stavR_dtstars.stas); stav_dtstars(ia) = stavR_dtstars.dtstar(ib); Nobs_dtstars(ia) = stavR_dtstars.Nobs(ib);

stav_ydT  = stav_dTp~=0 & stav_dTs~=0 & ~isnan(stav_dTp) & ~isnan(stav_dTs) & abs(stav_dTp)<3 & abs(stav_dTs)<5 & Nobs_dTp>4 & Nobs_dTs>4;
stav_ydts = stav_dtstarp~=0 & stav_dtstars~=0 & ~isnan(stav_dtstarp) & ~isnan(stav_dtstars) & abs(stav_dtstarp)<3 & abs(stav_dtstars)<5 & Nobs_dtstarp>4 & Nobs_dtstars>4;
stav_yP   = stav_dTp~=0 & ~isnan(stav_dTp) & stav_dtstarp~=0 & ~isnan(stav_dtstarp) & abs(stav_dTp)<3 & abs(stav_dtstarp)<3 & Nobs_dTp>4 & Nobs_dtstarp>4;
stav_yS   = stav_dTs~=0 & ~isnan(stav_dTs) & stav_dtstars~=0 & ~isnan(stav_dtstars) & abs(stav_dTs)<5 & abs(stav_dtstars)<5 & Nobs_dTs>4 & Nobs_dtstars>4;

% stav diff-tt line fit
tts_stav = stav_dTs(stav_ydT) - mean( stav_dTs(stav_ydT) );
ttp_stav = stav_dTp(stav_ydT) - mean( stav_dTp(stav_ydT) );
[ m_tt_stav] = fit_LSqEr( ttp_stav,tts_stav,1,1,1:0.01:5,[],10,1);
[ m_tto_stav,b_tto_stav] = fit_LSqEr( ttp_stav(stav_isob(stav_ydT)),tts_stav(stav_isob(stav_ydT)),0,1,1:0.01:5,-10:0.1:10,10,1);

% stav diff-tstar line fit
tss_stav = stav_dtstars(stav_ydts) - mean( stav_dtstars(stav_ydts) );
tsp_stav = stav_dtstarp(stav_ydts) - mean( stav_dtstarp(stav_ydts) );
[ m_ts_stav] = fit_LSqEr( tsp_stav,tss_stav,1,1,1:0.01:5,[],4,1);
[ m_tso_stav,b_tso_stav] = fit_LSqEr( tsp_stav(stav_isob(stav_ydts)),tss_stav(stav_isob(stav_ydts)),0,1,1:0.01:5,-10:0.1:10,4,1);

%% PLOT STAV DIFFERENTIAL TRAVEL TIME
figure(11), clf, set(gcf,'position',[200 200 510 580])
hold on

scatter(ttp_stav,tts_stav,50,'k','o','LineWidth',1.8)
scatter(ttp_stav(stav_isob(stav_ydT)),tts_stav(stav_isob(stav_ydT)),50,'b','o','LineWidth',1.8)

axis([-2.7 2.7 -3.7 3.7]); ax = axis;

plot(ax(1:2),m_tt_stav*ax(1:2),'r--','LineWidth',3)

text(-2.4,3.4,sprintf('\\textbf{Slope = %.2f}',m_tt_stav),'FontSize',18,'interpreter','latex')
text(-2.2,3.05,sprintf('\\textbf{$\\Rightarrow\\, \\delta \\ln{V_S}/\\delta \\ln{V_P} = %.2f$}',m_tt_stav/1.81),'FontSize',18,'interpreter','latex')

% obs slope
plot(ax(1:2), m_tto_stav*ax(1:2) + b_tto_stav,'b--','LineWidth',1)
text((ax(3)-b_tto_stav)/m_tto_stav+0.5,ax(3)+0.55,...
    sprintf('OBS Slope\\ = %.2f ',m_tto_stav),'FontSize',14,'interpreter','latex','color','b')
text((ax(3)-b_tto_stav)/m_tto_stav+0.8,ax(3)+0.28,...
    sprintf('$\\Rightarrow\\, \\delta \\ln{V_S}/\\delta \\ln{V_P} = %.2f$',m_tto_stav/1.81),'FontSize',14,'interpreter','latex','color','b')

set(gca,'fontsize',12,'color','none','LineWidth',2,'box','on')
% title('Ratio of $P$ and $S$ travel times residuals','FontSize',18,'interpreter','latex')
xlabel('Sta-Av $P$-wave $\delta T$ (s)','FontSize',20,'interpreter','latex')
ylabel('Sta-Av $S$-wave $\delta T$ (s)','FontSize',20,'interpreter','latex')

if ifsave
save2pdf(11,['dTS_vs_dTP_stav',obsstr],figdir); pause(0.1)
end
%% PLOT STAV DIFFERENTIAL T-STAR
figure(12), clf, set(gcf,'position',[300 200 510 580])
hold on

scatter(tsp_stav,tss_stav,50,'k','o','LineWidth',1.8)
scatter(tsp_stav(stav_isob(stav_ydts)),tss_stav(stav_isob(stav_ydts)),50,'b','o','LineWidth',1.8)

axis([-2.7 2.7 -3.7 3.7]); ax = axis;

plot(ax(1:2),m_ts_stav*ax(1:2),'r--','LineWidth',3)

text(-2.4,3.4,sprintf('\\textbf{Slope = %.2f}',m_ts_stav),'FontSize',18,'interpreter','latex')
text(-2.2,3.05,sprintf('\\textbf{$\\Rightarrow\\, Q_P/Q_S = %.2f$}',m_ts_stav/1.81),'FontSize',18,'interpreter','latex')

% obs slope
plot(ax(1:2), m_tso_stav*ax(1:2) + b_tso_stav,'b--','LineWidth',1)
text((ax(3)-b_tso_stav)/m_tso_stav+0.3,ax(3)+0.55,...
    sprintf('OBS Slope\\ = %.2f ',m_tso_stav),'FontSize',14,'interpreter','latex','color','b')
text((ax(3)-b_tso_stav)/m_tso_stav+0.6,ax(3)+0.28,...
    sprintf('$\\Rightarrow\\, Q_P/Q_S = %.2f$',m_tso_stav/1.81),'FontSize',14,'interpreter','latex','color','b')

set(gca,'fontsize',12,'color','none','LineWidth',2,'box','on')
% title('Ratio of $P$ and $S$ travel times residuals','FontSize',18,'interpreter','latex')
xlabel('Sta-Av $P$-wave $\Delta t^*$ (s)','FontSize',20,'interpreter','latex')
ylabel('Sta-Av $S$-wave $\Delta t^*$ (s)','FontSize',20,'interpreter','latex')

if ifsave
save2pdf(12,['dtstarS_vs_dtstarP_stav',obsstr],figdir); pause(0.1)
end

%% PLOT STAV S TT vs. T-STAR
figure(13), clf, set(gcf,'position',[300 200 500 500]), hold on
scatter(stav_dTs(stav_yS),stav_dtstars(stav_yS),40,'k','o','LineWidth',1.5)
scatter(stav_dTs(stav_isob&stav_yS),stav_dtstars(stav_isob&stav_yS),40,'b','o','LineWidth',1.5)

axis([-3.7 3.7 -3.7 3.7])
set(gca,'fontsize',12,'color','none','LineWidth',2,'box','on')
% title('Ratio of differential travel time to t-star for S','FontSize',18,'interpreter','latex')
xlabel('$\delta T_S$','FontSize',20,'interpreter','latex')
ylabel('$\Delta t^*_S$','FontSize',20,'interpreter','latex')

save2pdf(13,['dtstar_vs_dT_S_stav',obsstr],figdir); pause(0.1)

%% PLOT STAV S TT vs. T-STAR
figure(14), clf, set(gcf,'position',[300 200 500 500]), hold on
scatter(stav_dTp(stav_yP),stav_dtstarp(stav_yP),40,'k','o','LineWidth',1.5)
scatter(stav_dTp(stav_isob&stav_yP),stav_dtstarp(stav_isob&stav_yP),40,'b','o','LineWidth',1.5)

axis([-2 2 -1.7 1.7])
set(gca,'fontsize',12,'color','none','LineWidth',2,'box','on')
% title('Ratio of differential travel time to t-star for S','FontSize',18,'interpreter','latex')
xlabel('$\delta T_P$','FontSize',20,'interpreter','latex')
ylabel('$\Delta t^*_P$','FontSize',20,'interpreter','latex')

if ifsave
save2pdf(14,['dtstar_vs_dT_P_stav',obsstr],figdir); pause(0.1)
end

% %% =======================  OBS STATION AVERAGES ======================= %%
% [ stav_dTpo,~ ] = lsq_sta_evt( dTpo,0.01 );
% [ stav_dTso,~ ] = lsq_sta_evt( dTso,0.01 );
% [ stav_dtstarpo,~ ] = lsq_sta_evt( dtstpo,0.01 );
% [ stav_dtstarso,~ ] = lsq_sta_evt( dtstso,0.01 );
% 
% stav_ydTo  = stav_dTpo~=0   & stav_dTso~=0   & ~isnan(stav_dTpo)   & ~isnan(stav_dTso)   & abs(stav_dTpo)<3   & abs(stav_dTso)<5;
% stav_ydtso = stav_dtstarpo~=0 & stav_dtstarso~=0 & ~isnan(stav_dtstarpo) & ~isnan(stav_dtstarso) & abs(stav_dtstarpo)<3 & abs(stav_dtstarso)<5;
% stav_yPo   = stav_dTpo~=0 & ~isnan(stav_dTpo) & stav_dtstarpo~=0 & ~isnan(stav_dtstarpo) & abs(stav_dTpo)<3 & abs(stav_dtstarpo)<3;
% stav_ySo   = stav_dTso~=0 & ~isnan(stav_dTso) & stav_dtstarso~=0 & ~isnan(stav_dtstarso) & abs(stav_dTso)<5 & abs(stav_dtstarso)<5;
% 
% % stav diff-tt line fit
% ttso_stav = stav_dTso(stav_ydTo) - mean( stav_dTso(stav_ydTo) );
% ttpo_stav = stav_dTpo(stav_ydTo) - mean( stav_dTpo(stav_ydTo) );
% [ m_tto_stav] = fit_LSqEr( ttpo_stav,ttso_stav,1,1,1:0.01:5,[],10,1);
% 
% % stav diff-tstar line fit
% tsso_stav = stav_dtstarso(stav_ydtso) - mean( stav_dtstarso(stav_ydtso) );
% tspo_stav = stav_dtstarpo(stav_ydtso) - mean( stav_dtstarpo(stav_ydtso) );
% [ m_tso_stav] = fit_LSqEr( tspo_stav,tsso_stav,1,1,1:0.01:5,[],4,1);
% 
% %% PLOT OBS STAV DIFFERENTIAL TRAVEL TIME
% figure(15), clf, set(gcf,'position',[300 200 480 610])
% hold on
% 
% scatter(ttpo_stav,ttso_stav,40,'b','o','LineWidth',1.5)
% 
% axis([-3.7 3.7 -5 5]); ax = axis;
% 
% plot(ax(1:2),m_tto_stav*ax(1:2),'r--','LineWidth',3)
% 
% text(-3.2,4.6,sprintf('\\textbf{Slope = %.2f}',m_tto_stav),'FontSize',17,'interpreter','latex')
% text(-3,4.15,sprintf('\\textbf{$\\Rightarrow\\, \\delta \\ln{V_S}/\\delta \\ln{V_P} = %.2f$}',m_tto_stav/1.81),'FontSize',17,'interpreter','latex')
% 
% set(gca,'fontsize',12,'color','none','LineWidth',2,'box','on')
% % title('Ratio of $P$ and $S$ travel times residuals','FontSize',18,'interpreter','latex')
% xlabel('OBS Sta-Av $P$ residual time (s)','FontSize',16,'interpreter','latex')
% ylabel('OBS Sta-Av $S$ residual time (s)','FontSize',16,'interpreter','latex')
% 
% if ifsave
% save2pdf(15,'dTS_vs_dTP_stav_OBS',figdir); pause(0.1)
% end
% 
% %% PLOT OBS STAV DIFFERENTIAL T-STAR
% figure(16), clf, set(gcf,'position',[300 200 480 610])
% hold on
% 
% scatter(tspo_stav,tsso_stav,40,'b','o','LineWidth',1.5)
% 
% axis([-3.7 3.7 -5 5]); ax = axis;
% 
% plot(ax(1:2),m_tso_stav*ax(1:2),'r--','LineWidth',3)
% 
% text(-3.2,4.6,sprintf('\\textbf{Slope = %.2f}',m_tso_stav),'FontSize',17,'interpreter','latex')
% text(-3,4.15,sprintf('\\textbf{$\\Rightarrow\\, Q_P/Q_S = %.2f$}',m_tso_stav/1.81),'FontSize',17,'interpreter','latex')
% 
% set(gca,'fontsize',12,'color','none','LineWidth',2,'box','on')
% % title('Ratio of $P$ and $S$ travel times residuals','FontSize',18,'interpreter','latex')
% xlabel('OBS $P$-wave $\Delta t^*$ (s)','FontSize',16,'interpreter','latex')
% ylabel('OBS $S$-wave $\Delta t^*$ (s)','FontSize',16,'interpreter','latex')
% 
% if ifsave
% save2pdf(16,'dtstarS_vs_dtstarP_stav_OBS',figdir); pause(0.1)
% end
% 
% %% PLOT OBS STAV S TT vs. T-STAR
% figure(17)
% scatter(stav_dTso(stav_ySo),stav_dtstarso(stav_ySo),40,'b','o')
% 
% set(gca,'fontsize',12,'color','none','LineWidth',2,'box','on')
% % title('Ratio of differential travel time to t-star for S','FontSize',18,'interpreter','latex')
% xlabel('OBS $\delta T_S$','FontSize',15,'interpreter','latex')
% ylabel('OBS $\Delta t^*_S$','FontSize',15,'interpreter','latex')
% 
% %% PLOT OBS STAV S TT vs. T-STAR
% figure(18)
% scatter(stav_dTpo(stav_yPo),stav_dtstarpo(stav_yPo),40,'b','o')
% 
% set(gca,'fontsize',12,'color','none','LineWidth',2,'box','on')
% % title('Ratio of differential travel time to t-star for S','FontSize',18,'interpreter','latex')
% xlabel('OBS $\delta T_P$','FontSize',15,'interpreter','latex')
% ylabel('OBS $\Delta t^*_P$','FontSize',15,'interpreter','latex')

end % if stav

%% ======================  KEY TO ALL THE SAVING ======================= %%
% save2pdf(1 ,'dTs_vs_dTp_indiv',figdir); pause(0.1)
% save2pdf(11,'dTs_vs_dTp_stav',figdir); pause(0.1)
% save2pdf(5 ,'dTs_vs_dTp_indiv_OBS',figdir); pause(0.1)
% save2pdf(15,'dTs_vs_dTp_stav_OBS',figdir); pause(0.1)
% save2pdf(4 ,'dtstar_vs_dT_P_indiv',figdir); pause(0.1)
% save2pdf(14,'dtstar_vs_dT_P_stav',figdir); pause(0.1)
% save2pdf(3 ,'dtstar_vs_dT_S_indiv',figdir); pause(0.1)
% save2pdf(13,'dtstar_vs_dT_S_stav',figdir); pause(0.1)
% save2pdf(2 ,'dtstars_vs_dtstarp_indiv',figdir); pause(0.1)
% save2pdf(12,'dtstars_vs_dtstarp_stav',figdir); pause(0.1)
% save2pdf(6 ,'dtstars_vs_dtstarp_indiv_OBS',figdir); pause(0.1)
% save2pdf(16,'dtstars_vs_dtstarp_stav_OBS',figdir); pause(0.1)

