% This script compares the traces and spectra for just two stations for a
% given event. It plots all three methods for obtaining relative spectra
% (specR, comb, fAxP) as well as the associated best fits. It chooses the
% pair of stations on the basis of how good the fits are, aiming to
% showcase the methods working at their best.
clear all
% close all

addpath('~/Dropbox/MATLAB/AttenBody/matguts','~/Dropbox/MATLAB/AttenBody/');

% project details
dbname = 'EARdb';
phase = 'P';
comp  = 'Z';
methods = {'specR','comb','fAxP'}'; % must be column vector
methcol = [27,158,119; 247,95,2; 117,112,179]./255; % colorbrewer (take rows)

%% random
spd = 24*3600;

%% Get to work
dbdir = ['/Users/zeilon/Dropbox/Work/',dbname,'/'];
cd(dbdir);

% look at events available, and select orid
reseval_2_EVT_RESULTS(dbname,dbdir,{phase},{comp})
orid = input('select preferred orid: ');
run([dbdir,dbname,'_startup.m'])

%% load relevant files
evt = dir([data_eqar_dir,num2str(orid),'*']);
% get eqar
load([data_eqar_dir,evt.name,'/_EQAR_',phase,'_',comp])

% bring specR notation in line
% bring "dtstar" vs. "specR" notation into alignment
if ~isfield(eqar,'dtstar_specR') && isfield(eqar,'dtstar')
    [eqar.dtstar_specR] = deal(eqar.dtstar);
    [eqar.dtstarstd_specR] = deal(eqar.std_dtstar);
    eqar = rmfield(eqar,{'dtstar','std_dtstar'});
end

% get pairs
prs_c = dir([resdir,'PAIRSPECS/',num2str(orid),'_pairspecs_comb_',phase,comp,'.mat']);
prs_f = dir([resdir,'PAIRSPECS/',num2str(orid),'_pairspecs_fAxP_',phase,comp,'.mat']);
prs_c = load([resdir,'PAIRSPECS/',prs_c.name]);
prs_f = load([resdir,'PAIRSPECS/',prs_f.name]);

%% associate stations with pairs
% comb
prs_c.pairwise.st1 = cell(length(prs_c.pairwise.dtstar),1);
prs_c.pairwise.st2 = cell(length(prs_c.pairwise.dtstar),1);
ipr = 0;
for is1 = 1:length(prs_c.sts)
    for is2 = is1+1:length(prs_c.sts)
        ipr = ipr+1;
        prs_c.pairwise.st1{ipr} = prs_c.sts{is1};
        prs_c.pairwise.st2{ipr} = prs_c.sts{is2};
    end, clear is2;
end, clear is1;
% fAxP
prs_f.pairwise.st1 = cell(length(prs_f.pairwise.dtstar),1);
prs_f.pairwise.st2 = cell(length(prs_f.pairwise.dtstar),1);
ipr = 0;
for is1 = 1:length(prs_f.sts)
    for is2 = is1+1:length(prs_f.sts)
        ipr = ipr+1;
        prs_f.pairwise.st1{ipr} = prs_f.sts{is1};
        prs_f.pairwise.st2{ipr} = prs_f.sts{is2};
    end, clear is2;
end, clear is1;

%% associate misfits
Nsta = length(eqar);
misfits = nan(Nsta,Nsta,2);
% comb
for is1 = 1:Nsta
    for is2 = is1+1:Nsta
        sta1 = eqar(is1).sta;
        sta2 = eqar(is2).sta;
        % comb pair with these two stas
        ipr_c = find(...
                    (strcmp(prs_c.pairwise.st1,sta1) & strcmp(prs_c.pairwise.st2,sta2)) |...
                    (strcmp(prs_c.pairwise.st1,sta2) & strcmp(prs_c.pairwise.st2,sta1)) );
        % fAxP pair with these two stas
        ipr_f = find(...
                    (strcmp(prs_f.pairwise.st1,sta1) & strcmp(prs_f.pairwise.st2,sta2)) |...
                    (strcmp(prs_f.pairwise.st1,sta2) & strcmp(prs_f.pairwise.st2,sta1)) );
        if ~isempty(ipr_c)
            misfits(is1,is2,1) = prs_c.pairwise.misfit_normed(ipr_c);
        end
        if ~isempty(ipr_f)
            misfits(is1,is2,2) = prs_f.pairwise.misfit_normed(ipr_f);
        end
    end, clear is2;
end, clear is1;

% find the preferred point (only care about fAxP, inclue 1e-9 comb to
% ensure not nan)
[~,is2_best,is1_best] = mingrid(misfits(:,:,2) + 1e-9*misfits(:,:,2));

eqar1 = eqar(is1_best); % just subset to layer for sta1
eqar2 = eqar(is2_best); % just subset to layer for sta2

sta1 = eqar1.sta;
sta2 = eqar2.sta;

% comb pair with these two stas
ipr_c = find(...
            (strcmp(prs_c.pairwise.st1,sta1) & strcmp(prs_c.pairwise.st2,sta2)) |...
            (strcmp(prs_c.pairwise.st1,sta2) & strcmp(prs_c.pairwise.st2,sta1)) );
% fAxP pair with these two stas
ipr_f = find(...
            (strcmp(prs_f.pairwise.st1,sta1) & strcmp(prs_f.pairwise.st2,sta2)) |...
            (strcmp(prs_f.pairwise.st1,sta2) & strcmp(prs_f.pairwise.st2,sta1)) );



%% NOW move on with these two stations
    
%% ========== Load the two traces ==========
dat1 = eqar1.(['dat',comp]);
dat2 = eqar2.(['dat',comp]);

tt1 = (eqar1.tt'- eqar1.abs_arrT)*spd;% - eqar1.dT;
tt2 = (eqar2.tt'- eqar2.abs_arrT)*spd;% - eqar2.dT;

% traces for pre-shifted dat2
tt2_ = (eqar2.tt'- eqar2.abs_arrT)*spd + (eqar2.dT - eqar1.dT);
dat2_ = dat2;  dat2_(abs(tt2_)>5)= nan;

%% predict differential spectra
frq = eqar1.frq;
ws = 2*pi*frq;

%% make structure with info for all three methods
gdf_c = prs_c.pairwise.wts(ipr_c,:)~=0;
gdf_f = prs_f.pairwise.wts(ipr_f,:)~=0;

SCF = struct('method',methods...
            ,'freq',...
                    {eqar1.frq;prs_c.fmids(gdf_c);prs_f.fmids(gdf_f)}...
            ,'lnA',...
                    {log(eqar2.specs./eqar1.specs);
                     log(prs_c.pairwise.As(ipr_c,gdf_c))';
                     log(prs_f.pairwise.As(ipr_f,gdf_f))'}...
            ,'phi',...
                    {nan;
                     (prs_c.pairwise.phis(ipr_c,gdf_c))';
                     (prs_f.pairwise.phis(ipr_f,gdf_f))'}...
            ,'wts',...
                    {ones(size(eqar1.frq));
                     (prs_c.pairwise.wts(ipr_c,gdf_c))';
                     (prs_f.pairwise.wts(ipr_f,gdf_f))'}...
            ,'prex',...
                    {eqar1.par_dtstar_specR.swindow(1);
                     -eqar1.par_dtstar_comb.wind.prex;
                     -eqar1.par_dtstar_fAxP.wind.prex}...
            ,'postx',...
                    {eqar1.par_dtstar_specR.swindow(2);
                     eqar1.par_dtstar_comb.wind.postx;
                     eqar1.par_dtstar_fAxP.wind.postx}...
            ,'alpha',...
                    {0;eqar1.alpha_comb;eqar1.alpha_fAxP}...
            ,'amp2phiwt',...
                    {nan;
                     eqar1.par_dtstar_comb.inv.amp2phiwt;
                     eqar1.par_dtstar_fAxP.inv.amp2phiwt}...
            );
% dtstar
for im = 1:length(methods), imethod = methods{im};
    SCF(im).dtstar = eqar2.(['dtstar_',imethod]) - eqar1.(['dtstar_',imethod]);
    SCF(im).dtstarstd = eqar2.(['dtstarstd_',imethod]) + eqar1.(['dtstarstd_',imethod]);
end
% A0 and dT
SCF(1).dA0 = eqar2.specss(1)./eqar1.specss(1);
SCF(1).dT = nan;
for im = 2:3, imethod = methods{im};
    SCF(im).dA0 = eqar2.(['A0_',imethod])./eqar1.(['A0_',imethod]);
    SCF(im).dT = eqar2.(['dT_',imethod]) - eqar1.(['dT_',imethod]);
end


%% =========================================================================
%% ====================== START THE PLOTS ==================================

figure(130), clf, set(gcf,'pos',[600 309 736 1036])
ax1 = axes('pos',[0.1100 0.8 0.800 0.17]); hold on
ax2 = axes('pos',[0.1100 0.42 0.800 0.31]); hold on
ax3 = axes('pos',[0.1100 0.08 0.800 0.31]); hold on

% ================ PLOT THE WAVEFORMS ================
plot(ax1,tt1,dat1./max(abs(dat1)),'k','LineWidth',2)
plot(ax1,tt2,dat2./max(abs(dat1)),'r','LineWidth',2)
plot(ax1,tt2_,dat2_./max(abs(dat2)),'r:','LineWidth',2)

% dT window
% plot(ax1,-eqar1.par_dT.prex*[1 1],[-1 1],'b--',eqar1.par_dT.postx*[1 1],[-1 1],'b--','linewidth',1.5)
for im = 1:length(methods)
    hwind(im) = patch(ax1,[SCF(im).prex,SCF(im).prex,SCF(im).postx,SCF(im).postx,SCF(im).prex]',...
              0.97 + 0.08*im + 0.035*[1 -1 -1 1 1]',...
              methcol(im,:),'facealpha',0.6,'linestyle','none');
end


% station names, annotations
text(ax1,min([SCF.prex])+5,-0.5,sprintf('\\textbf{%s}',sta1),'interpreter','latex','fontsize',12,'horizontalalignment','left','verticalalignment','bottom','color','k')
text(ax1,min([SCF.prex])+5,-0.8,sprintf('\\textbf{%s}',sta2),'interpreter','latex','fontsize',12,'horizontalalignment','left','verticalalignment','bottom','color','r')
text(ax1,min([SCF.prex])-5,-1.65,'Time $\rightarrow$','interpreter','latex','fontsize',15,'horizontalalignment','left','verticalalignment','bottom')


% xlabel('Time, s','FontSize',22,'interpreter','latex')
set(ax1,'xlim',[min([SCF.prex])-10 max([SCF.postx])+15],'ylim',[-1.1 1.28],...%'XAxisLocation','top',...
    'Fontsize',12,'linewidth',2,'box','on','ytick',[])



% ================ PLOT THE AMPLITUDE SPECTRA ================
% spectral ratio measurements
plot(ax2,SCF(1).freq,SCF(1).lnA,':s','linewidth',1.5,'color',methcol(1,:))
% comb and fAxP measurements
for im = 2:length(methods)
    ha(im) = scatter(ax2,SCF(im).freq,SCF(im).lnA,100*SCF(im).wts,...
        methcol(im,:),'filled','markeredgecolor','k','linewidth',1.5,'markerfacealpha',0.5,'MarkerEdgeAlpha',0.5);
end

% ================ PLOT THE PHASE SPECTRA ================
% comb and fAxP measurements
for im = 2:length(methods)
    hp(im) = scatter(ax3,SCF(im).freq,SCF(im).phi,100*SCF(im).wts,...
        methcol(im,:),'filled','markeredgecolor','k','linewidth',1.5,'markerfacealpha',0.5,'MarkerEdgeAlpha',0.5);
end


% predictions - best fits from all station averages
for im = 2:3

    % predict spectra from event-wide (eqar) value
    if SCF(im).alpha == 0
        lnApred = log(SCF(im).dA0) - pi*SCF(im).freq*SCF(im).dtstar;
        phipred = -(1/pi)*log(SCF(im).freq)*SCF(im).dtstar + SCF(im).dtstar + SCF(im).dT;
    else
        w0 = 2*pi;
        lnApred = log(SCF(im).dA0) - 0.5 * ws.^(1-SCF(im).alpha) * w0^SCF(im).alpha * SCF(im).dtstar;
        phipred = 0.5*ws.^(-SCF(im).alpha).*(2*pi).^SCF(im).alpha * cot(SCF(im).alpha*pi/2)*SCF(im).dtstar + SCF(im).dT;
    end    

    % amp prediction
    hap(im) = plot(ax2,frq,lnApred,...
        'color',methcol(im,:),'linewidth',1.5);
    % phi prediction
    hpp(im) = plot(ax3,frq,phipred,...
        'color',methcol(im,:),'linewidth',1.5);

    % fits and predict spectra from this station pair alone
    [ dtstar_bf(im),dT_bf(im),A0_bf(im)] ...
        = invert_1pair_Aphi_4_dtdtstar( exp(SCF(im).lnA),SCF(im).phi,SCF(im).freq,SCF(im).wts,SCF(im).amp2phiwt,SCF(im).alpha);

    lnAfit = log(A0_bf(im)) - 0.5 * ws.^(1-SCF(im).alpha) * w0^SCF(im).alpha * dtstar_bf(im);
    phifit = 0.5*ws.^(-SCF(im).alpha).*(2*pi).^SCF(im).alpha * cot(SCF(im).alpha*pi/2)*dtstar_bf(im) + dT_bf(im);
    
    % amp fit
    haf(im) = plot(ax2,frq,lnAfit,'--',...
        'color',methcol(im,:),'linewidth',1.5);
    % phi fit
    hpf(im) = plot(ax3,frq,phifit,'--',...
        'color',methcol(im,:),'linewidth',1.5);
    
end

% accentuate fAxP plots
set([ha(3),hp(3)],'MarkerFaceAlpha',0.8);
set([hap(3),hpp(3)],'LineWidth',2.5)


%
% specR
% for imethod = 1:3
%         [ dtstar,dT,A0,misfit,resnorm,A_SSE,P_SSE,lnA_SSE,dtstar_std,dT_std,dlnA0_std ] ...
%         = invert_1pair_Aphi_4_dtdtstar( As(inds),phis(inds),ff_use(inds), wts(inds),amp2phiwt,alpha);
% 
%     for alp = [0,eqar1.par_dtstar_comb.inv.alpha]
%         if alp == 0
%             lnApred = log(A0) - pi*fmids*dtstar;
%             phipred = -(1/pi)*log(fmids)*dtstar + dtstar + dT;
%         else
%             lnApred = log(A0) - 0.5 * ws.^(1-alp) * w0^alp * dtstar;
%             phipred = 0.5*ws.^(-alp).*(2*pi).^alp * cot(alp*pi/2)*dtstar + dT;
%         end
% 
%     end
% end

% 
% plot(fmids,lnApred,'g','Linewidth',2.5)
% plot(fmids,lnApred_a27,'--g','Linewidth',1.5)
% % plot predictions from all evts...
% plot(fmids,lnApred_specR,'r','Linewidth',1.5)
% plot(fmids,lnApred_comb,'c','Linewidth',1.5)
% 
% plot(ff_fft,specRfft ,'-b','linewidth',2)
% scatter(ff_fft,specRfft,40*wt_fft,'markerfacecolor','b')
% plot(ff_fft,lnApred_fft,'--b','Linewidth',1.5)
% 
% gdf = wts~=0;
% scatter(fmids(gdf),log(As(gdf)),150*wts(gdf),'or','MarkerFaceColor','r')
% scatter(fmids(inds),log(As(inds)),150*wts(inds),'ok','MarkerFaceColor','r','linewidth',1.5)

% ha1 = plot([eqar(ind1).fcross,1],[eqar(ind1).fcross,0],1.9,90,4,0.12,'FaceColor','k'); % fcross for sta 1
% ha2 = arrow([eqar(ind2).fcross,1],[eqar(ind2).fcross,0],1.9,90,4,0.12,'FaceColor','k'); % fcross for sta 2
% text(eqar(ind1).fcross,1.3,'$f^{max}$','interpreter','latex','fontsize',18,'horizontalalignment','center','verticalalignment','bottom')
% text(eqar(ind2).fcross,1.3,'$f^{max}_2$','interpreter','latex','fontsize',16,'horizontalalignment','center','verticalalignment','bottom')
xline(ax2,min([eqar1.fcrosshi,eqar2.fcrosshi]),'--b','LineWidth',1.5)
xline(ax2,max([eqar1.fcrosslo,eqar2.fcrosslo]),'--b','LineWidth',1.5)

% xlabel('freq (Hz)','FontSize',16,'interpreter','latex')


xlabel(ax3,'frequency (Hz)','FontSize',22,'interpreter','latex')
ylabel(ax3,'$\Delta \psi_{12}$ (s)','FontSize',22,'interpreter','latex')
ylabel(ax2,'$\ln\,(R_{12})$','FontSize',22,'interpreter','latex')

set(ax2,'ylim',[-3.5 -0.5],'fontsize',15,'Xscale','linear','xticklabel',[],'linewidth',2,'box','on')
set(ax3,'fontsize',15,'Xscale','linear','linewidth',2,'box','on')
xlim([ax2,ax3],[0.00 1.2])

return

% ================ PLOT THE PHASE SPECTRA ================
% subplot(5,1,4:5), hold on
% plot(fmids,phipred,'g','Linewidth',2.5)
% plot(fmids,phipred_a27,'--g','Linewidth',1.5)
% % plot predictions from all evts...
% plot(fmids,phipred_specR,'r','Linewidth',1.5)
% plot(fmids,phipred_comb,'c','Linewidth',1.5)
% 
% plot(ff_fft,specPfft ,'-b','linewidth',2)
% scatter(ff_fft,specPfft,40*wt_fft,'markerfacecolor','b')
% plot(ff_fft,phipred_fft ,'--b','linewidth',2)

% scatter(fmids(gdf),phis(gdf),150*wts(gdf),'or','MarkerFaceColor','r')
% scatter(fmids(inds),phis(inds),150*wts(inds),'ok','MarkerFaceColor','r','linewidth',1.5)

% ha1 = arrow([eqar(ind1).fcross,1],[eqar(ind1).fcross,0],5,90,2,0.1,'FaceColor','m'); % fcross for sta 1
% ha2 = arrow([eqar(ind2).fcross,1],[eqar(ind2).fcross,0],5,90,2,0.1,'FaceColor','g'); % fcross for sta 2
% plot(fmax*[1 1],[-1 1],'--b')  
text(0.5,-.2,['$\mathbf{\Delta t^* = ',num2str(dtstar,'%.1f'),'}$'],'interpreter','latex','fontsize',20,'horizontalalignment','right','verticalalignment','bottom')

text(0.03,0,'More delayed $\rightarrow$','interpreter','latex','fontsize',16,...
    'rotation',90,'horizontalalignment','left','verticalalignment','bottom')



% ================ SAVE FIGURE ================
if ifsave
    ofile = sprintf('eg_comb_orid%.0f_%s%s_%s_v_%s_alp%.2f',orid,phase,comp,sta1,sta2,alp);
    save2pdf(30,ofile,'figs');
end



% samprate = eqar(ind1).samprate;
% dt = 1./samprate;
% fnq = samprate/2;
% T = prex+postx;
% fmax = mean([eqar([ind1,ind2]).fcross]);
% 
% %interp to twin
% tt = [-pretime:dt:pretime-dt]';
% dat1 = interp1(tt1,dat1,tt);
% dat2 = interp1(tt2,dat2,tt);
% 
% %% MAKE STACK INSTEAD
% if ischar(ind2)
% indgd = find(~isnan([eqar.dT]));
% tt2 = [-pretime:dt:pretime]'; tt2 = tt2(1:length(tt1));
% for ig = 1:length(indgd)
%     is = indgd(ig);
%     att = eqar(is).tt-eqar(is).abs_arrT; % shift to since MEASURED absolute arrival
%     all_dat0(:,is) = interp1(att,eqar(is).datT,tt2,'linear',0);
% end % loop on stas
% indgd = 1:size(eqar);
% indgd(mean(abs(all_dat0(:,indgd)))<1e-6)     = []; % kill zero traces
% indgd(isnan(mean(abs(all_dat0(:,indgd))))) = []; % kill nan traces
% stak = sum(all_dat0(:,indgd),2)/length(indgd);
% dat2 = stak;
% end
% 
% [ dat1 ] = filt_quick( dat1,1./40,10,dt);
% [ dat2 ] = filt_quick( dat2,1./40,10,dt);
% 
% %% Plot two traces
% figure(1), clf, set(gcf,'pos',[100 550 1000 350]), hold on
% plot(tt,dat1./max(abs(dat1)),'k','LineWidth',1.5)
% plot(tt,dat2./max(abs(dat2)),'r','LineWidth',1.5)
% plot(-prex*[1 1],max(abs(dat1))*[-1 1],'b--',postx*[1 1],max(abs(dat1))*[-1 1],'b--')
% 
% %% Make set of period windows for bandpass filter
% Tmids = logspace(log10(Tmin),log10(Tmax),Nwds)';
% Twdhs = 0.5*diff(logspace(log10(Tmin/2),log10(2*Tmax),Nwds+1)');
% % Twdhs = 0.5*linspace(1,3,Nwds);
% fmids = 1./Tmids;
% 
% 
% 
% %% Test combspectra function
% parms.comb.Tmin = Tmin; % fnq/2
% parms.comb.Tmax = Tmax;
% parms.comb.Nwds = Nwds;
% parms.comb.Tw_opt = 'scale';
% parms.comb.npol = 4;
% 
% parms.wind.pretime = pretime;
% parms.wind.prex = prex;
% parms.wind.postx = postx;
% parms.wind.taperx = 0.1;
% 
% parms.qc.minacor = minacor;
% parms.qc.maxphi = maxphi;
% parms.qc.mtmRcomp = true; % option to compare comb and mtm estimates of As
% parms.qc.focus_threshold = 3; % threshold value for focus metric to flag up as bad
% 
% parms.inv.amp2phiwt = amp2phiwt;
% parms.inv.fmin = .1;
% parms.inv.fmax = fmax;
% parms.inv.ifwt = true;
% parms.inv.corr_c_skip = true;
% parms.inv.fmin_cskip = 0.1;
% parms.inv.R2default = 0.9; % if zero, do not weight by R2 of each pairwise fit, Else, weight by multiple of this default.
% parms.inv.opt = 1; %USE 1!   for calc_fdependent...  1 is all in one method; 2 is one-by-one method
% 
% parms.inv.alpha = alp;
% 
% % GET station-by-station acceptable f range from crossing freqs
% fcross = [[eqar([ind1,ind2]).fcrosslo]',[eqar([ind1,ind2]).fcrosshi]'];
% fcross(:,2) = 0.3;
% ifplot = false;
% % 
% % tic
% % [delta_tstar,delta_T,std_dtstar,pairwise] = combspectra([dat1,dat2],samprate,parms,1);
% % fprintf('combspectra dtstar est\t= %.3f \n',diff(delta_tstar))
% % fprintf('combspectra dT est\t= %.3f \n',diff(delta_T))
% % toc
% % tic
% % [delta_tstar_skip,delta_T_skip,std_dtstar_skip,pairwise] = combspectra([dat1,dat2],fcross,samprate,eqar(1).par_dtstar_comb,1);
% [delta_tstar_skip,delta_T_skip,std_dtstar_skip,pairwise,fmids] = combspectra([dat1,dat2],fcross,samprate,parms,combspectraplotopt);
% dtstar_specR = eqar(ind2).dtstar - eqar(ind1).dtstar;
% dtstar_comb  = eqar(ind2).dtstar_comb - eqar(ind1).dtstar_comb;
% 
% fprintf('whole event dtstar_specR est\t= %.3f \n',dtstar_specR)
% fprintf('whole event dtstar_comb est\t= %.3f \n',dtstar_comb)
% fprintf('whole event dT est    \t= %.3f \n',eqar(ind2).dT-eqar(ind1).dT)
% fprintf('BETTER? combspectra dtstar est\t= %.3f \n',diff(delta_tstar_skip))
% fprintf('BETTER? combspectra dT est\t= %.3f \n',diff(delta_T_skip))
% % toc
% 
% [gdofgd,focusscore] = amp2phi_QC_focusing(parms,pairwise,fmids,1,{sta1,sta2});
% fprintf('Focus score:\n')
% disp(focusscore)
% 
% %% Get ready to plot
% As = pairwise.As';
% phis = pairwise.phis';
% wts = pairwise.wts';
% inds = find(fmids<=fmax & abs(phis)<maxphi & sqrt(wts)>minacor);
% [ dtstar,dT,A0,misfit,res ] = invert_1pair_Aphi_4_dtdtstar( As(inds),phis(inds),fmids(inds), wts(inds),amp2phiwt,alp);
% frq = eqar(ind1).frq;
% 
% %% NEW - diff. phase and amplitude spectra from fourier decomp. 
% % subset data to signal windows
% cp = struct('samprate',1./dt,'pretime',pretime,'prex',prex,'postx',postx,...
%                 'taperx',0.1,'fhi',100,'flo',0,'npoles',2,'norm',0);
% [ dat1w ] = data_clean( dat1,cp); 
% [ dat2w ] = data_clean( dat2,cp);
% 
% [dat1_fft,ff_fft] = fft_ze(dat1w,dt);
% [dat2_fft,~]   = fft_ze(dat2w,dt);
% N = length(dat1_fft);
% fplt = (ff_fft>=max(fcross(:,1))) & (ff_fft<=min(fcross(:,2)));
% specRfft = smooth(log(abs(dat2_fft)./abs(dat1_fft)));
% % specPfft = smooth(unwrap(angle(dat2_fft./dat1_fft))./(-2*pi*ff_fft));
% specPfft = smooth((angle(dat2_fft./dat1_fft))./(-2*pi*ff_fft));
% 
% specRfft = specRfft(fplt);
% specPfft = specPfft(fplt);
% ff_fft = ff_fft(fplt);
% % wt_fft = 1./(0.5+abs(interp1(eqar(ind1).frq,log(eqar(ind2).specs./eqar(ind1).specs),ff_fft) - specRfft));
% wt_fft_a = 0.4./(.2+abs(interp1(eqar(ind1).frq,log(eqar(ind2).specs./eqar(ind1).specs),ff_fft) - specRfft));
% wt_fft_p = [1;abs(diff(specPfft)./diff(ff_fft))<10];
% wt_fft = (wt_fft_a.*wt_fft_p) + 0.0001;
% figure(444),clf
% subplot(211), hold on
% plot(ff_fft,specRfft ,'-b','linewidth',2)
% plot(eqar(ind1).frq,log(eqar(ind2).specs./eqar(ind1).specs),'-r')
% plot(ff_fft,wt_fft,'k')
% xlim([0 0.5])
% subplot(212), hold on
% plot(ff_fft,specPfft ,'-b','linewidth',2)
% xlim([0 0.5])
% 
% 
% 
% 
% %% best-fit f dependence
% alpha_test = [0:0.05:0.9]';
% for ia = 1:length(alpha_test)
%     [ alpha_dtstar(ia),~,~,alpha_misfit(ia),~ ] = invert_1pair_Aphi_4_dtdtstar( As(inds),phis(inds),fmids(inds), wts(inds),amp2phiwt,alpha_test(ia));
% end
% figure(11), clf
% plot(alpha_test,alpha_misfit,'o-')
% figure(12), clf
% plot(alpha_test,alpha_dtstar,'o-')
% 
% %% predicted slopes from this pair of stas
% ws = 2*pi*fmids;
% if alp == 0
%     lnApred   = log(A0) - pi*fmids*dtstar;
%     phipred = -(1/pi)*log(fmids)*dtstar + dtstar + dT;
% else
%     lnApred = log(A0) - 0.5 * ws.^(1-alp) * w0^alp * dtstar;
%     phipred = 0.5*ws.^(-alp).*(2*pi).^alp * cot(alp*pi/2)*dtstar + dT;
% end
% 
% %% predicted slopes from inversion of all stas - in eqar
% % specR
% lnApred_specR   = log(A0) - pi*fmids*dtstar_specR;
% phipred_specR = -(1/pi)*log(fmids)*dtstar_specR + dtstar_specR + dT;
% 
% % comb
% lnApred_comb   = log(A0) - pi*fmids*dtstar_comb;
% phipred_comb = -(1/pi)*log(fmids)*dtstar_comb + dtstar_comb + dT;
% 
% 
% %% do for alp0.27 too
% [ dtstar_a27,dT_a27,A0_a27,misfit,res ] = invert_1pair_Aphi_4_dtdtstar( As(inds),phis(inds),fmids(inds), wts(inds),amp2phiwt,0.27);
% lnApred_a27 = log(A0_a27) - 0.5 * ws.^(1-0.27) * w0^0.27 * dtstar_a27;
% phipred_a27 = 0.5*ws.^(-0.27).*(2*pi).^0.27 * cot(0.27*pi/2)*dtstar_a27 + dT_a27;
% 
% %% do for fft spectra too
% [ dtstar_fft,dT_fft,A0_fft] = invert_1pair_Aphi_4_dtdtstar( exp(specRfft),specPfft,ff_fft, wt_fft,amp2phiwt,0);
% lnApred_fft = log(A0_fft) - pi*ff_fft * dtstar_fft;
% phipred_fft = -(1/pi)*log(ff_fft)*dtstar_fft + dtstar_fft + dT_fft;
% fprintf('BETTER? specXPR dtstar est = %.3f \n',dtstar_fft)
% fprintf('BETTER? specXPR dT est\t   = %.3f \n',dT_fft)
% 
% 
%     
% %% =========================================================================
% %% ====================== START THE PLOTS ==================================
% 
% figure(30), clf, set(gcf,'pos',[600 309 736 1036])
% 
% % ================ PLOT THE WAVEFORMS ================
% subplot(5,1,1), hold on
% plot(tt,dat1./max(abs(dat1)),'k','LineWidth',2)
% plot(tt,dat2./max(abs(dat2)),'r','LineWidth',2)
% plot(-prex*[1 1],[-1 1],'b--',postx*[1 1],[-1 1],'b--','linewidth',1.5)
% text(-prex-18,-0.5,sprintf('\\textbf{%s}',sta1),'interpreter','latex','fontsize',12,'horizontalalignment','left','verticalalignment','bottom','color','k')
% text(-prex-18,-0.8,sprintf('\\textbf{%s}',sta2),'interpreter','latex','fontsize',12,'horizontalalignment','left','verticalalignment','bottom','color','r')
% text(postx+15,-1.65,'Time $\rightarrow$','interpreter','latex','fontsize',15,'horizontalalignment','left','verticalalignment','bottom')
% 
% 
% % xlabel('Time, s','FontSize',22,'interpreter','latex')
% set(gca,'xlim',[-prex-20 postx+20],'ylim',[-1.1 1.1],...%'XAxisLocation','top',...
%     'Fontsize',12,'linewidth',2,'box','on','ytick',[])
% 
% % ================ PLOT THE AMPLITUDE SPECTRA ================
% subplot(5,1,2:3), hold on
% % plot(frq(2:end),log(spec2./spec1),'-ok')
% plot(eqar(ind1).frq,log(eqar(ind2).specs./eqar(ind1).specs),'-xk','linewidth',1.5)
% % plot(eqar(ind1).frq,log(eqar(ind2).specs./eqar(ind2).specn),'-ob')
% % plot(eqar(ind1).frq,log(eqar(ind1).specs./eqar(ind1).specn),'-om')
% plot(fmids,lnApred,'g','Linewidth',2.5)
% plot(fmids,lnApred_a27,'--g','Linewidth',1.5)
% % plot predictions from all evts...
% plot(fmids,lnApred_specR,'r','Linewidth',1.5)
% plot(fmids,lnApred_comb,'c','Linewidth',1.5)
% 
% plot(ff_fft,specRfft ,'-b','linewidth',2)
% scatter(ff_fft,specRfft,40*wt_fft,'markerfacecolor','b')
% plot(ff_fft,lnApred_fft,'--b','Linewidth',1.5)
% 
% gdf = wts~=0;
% scatter(fmids(gdf),log(As(gdf)),150*wts(gdf),'or','MarkerFaceColor','r')
% scatter(fmids(inds),log(As(inds)),150*wts(inds),'ok','MarkerFaceColor','r','linewidth',1.5)
% 
% % ha1 = plot([eqar(ind1).fcross,1],[eqar(ind1).fcross,0],1.9,90,4,0.12,'FaceColor','k'); % fcross for sta 1
% % ha2 = arrow([eqar(ind2).fcross,1],[eqar(ind2).fcross,0],1.9,90,4,0.12,'FaceColor','k'); % fcross for sta 2
% % text(eqar(ind1).fcross,1.3,'$f^{max}$','interpreter','latex','fontsize',18,'horizontalalignment','center','verticalalignment','bottom')
% % text(eqar(ind2).fcross,1.3,'$f^{max}_2$','interpreter','latex','fontsize',16,'horizontalalignment','center','verticalalignment','bottom')
% plot(max(fcross(:,1))*[1 1],[-1 1],'--b','LineWidth',1.5)
% plot(min(fcross(:,2))*[1 1],[-1 1],'--b','LineWidth',1.5)
% 
% % xlabel('freq (Hz)','FontSize',16,'interpreter','latex')
% ylabel('$\ln\,(R_{12})$','FontSize',22,'interpreter','latex')
% set(gca,'fontsize',15,'Xscale','linear','xticklabel',[],'xlim',[0.00 0.6],'ylim',[-3.5 2.5],'linewidth',2,'box','on')
% 
% % ================ PLOT THE PHASE SPECTRA ================
% subplot(5,1,4:5), hold on
% plot(fmids,phipred,'g','Linewidth',2.5)
% plot(fmids,phipred_a27,'--g','Linewidth',1.5)
% % plot predictions from all evts...
% plot(fmids,phipred_specR,'r','Linewidth',1.5)
% plot(fmids,phipred_comb,'c','Linewidth',1.5)
% 
% plot(ff_fft,specPfft ,'-b','linewidth',2)
% scatter(ff_fft,specPfft,40*wt_fft,'markerfacecolor','b')
% plot(ff_fft,phipred_fft ,'--b','linewidth',2)
% 
% scatter(fmids(gdf),phis(gdf),150*wts(gdf),'or','MarkerFaceColor','r')
% scatter(fmids(inds),phis(inds),150*wts(inds),'ok','MarkerFaceColor','r','linewidth',1.5)
% 
% % ha1 = arrow([eqar(ind1).fcross,1],[eqar(ind1).fcross,0],5,90,2,0.1,'FaceColor','m'); % fcross for sta 1
% % ha2 = arrow([eqar(ind2).fcross,1],[eqar(ind2).fcross,0],5,90,2,0.1,'FaceColor','g'); % fcross for sta 2
% % plot(fmax*[1 1],[-1 1],'--b')  
% text(0.5,-.2,['$\mathbf{\Delta t^* = ',num2str(dtstar,'%.1f'),'}$'],'interpreter','latex','fontsize',20,'horizontalalignment','right','verticalalignment','bottom')
% 
% text(0.03,0,'More delayed $\rightarrow$','interpreter','latex','fontsize',16,...
%     'rotation',90,'horizontalalignment','left','verticalalignment','bottom')
% 
% 
% xlabel('frequency (Hz)','FontSize',22,'interpreter','latex')
% ylabel('$\Delta \psi_{12}$ (s)','FontSize',22,'interpreter','latex')
% set(gca,'fontsize',15,'Xscale','linear','xlim',[0.00 0.6],'linewidth',2,'box','on')
% 
% 
% % ================ SAVE FIGURE ================
% if ifsave
%     ofile = sprintf('eg_comb_orid%.0f_%s%s_%s_v_%s_alp%.2f',orid,phase,comp,sta1,sta2,alp);
%     save2pdf(30,ofile,'figs');
% end