% Third script to test the method of measuring differential
% attenuation using combs of filters and measuring amplitude loss as well
% as phase lag within each narrow bandpass.
% This script uses real data and three different traces that have passed
% through different Q and V structure. 
% Second simple script to test the method of measuring differential
% attenuation using combs of filters and measuring amplitude loss as well
% as phase lag within each narrow bandpass.
% This script is attempts to do this with two traces that have passed
% through different Q and V material, to see how the velocity and Q
% differences trade off, and how the delta-tstar works.
clear all
close all
% close all
addpath('matguts')

datdir = '/Volumes/data/EastAfricaBW/';
phase = 'S';
comp  = 'T';


ifsave = false;

% filter combs parms
Tmin = 1;
Tmax = 25;
% Twid = 2
Nwds = 20;

% window parms
pretime = 100;
prex = 5;
postx = 20;

minacor = 0.2;
maxphi = 5;
amp2phiwt = 3;

alp = 0;
w0 = 2*pi;

fmax = 1;

combspectraplotopt = 1;

spd = 3600*24;



orid = 2200; % seaz ~ 299   
sta2 = 'E64'; % should be MORE attenuated
sta2 = 'AMME'; % should be MORE attenuated
% sta2 = 'DZEE'; % should be MORE attenuated
sta1 = 'E57'; % should be LESS attenuated
% sta1 = 'E67'; % should be LESS attenuated
% sta1 = 'E78'; % should be LESS attenuated


% orid = 819; % seaz ~ 299   
% sta2 = 'HALE'; % should be MORE attenuated
% sta1 = 'AHME'; % should be LESS attenuated

% orid = 881; % seaz ~ 299   
% sta2 = 'HITE'; % should be MORE attenuated
% sta1 = 'GALE'; % should be LESS attenuated

% orid = 1108; % seaz ~ 299   
% sta2 = 'MAYE'; % should be MORE attenuated
% sta1 = 'SAHE'; % should be LESS attenuated

% orid = 1094; % seaz ~ 299   
% sta2 = 'FINE'; % should be MORE attenuated
% sta1 = 'MAYE'; % should be LESS attenuated

% orid = 215; % seaz ~ 299  
% orid = 238; % seaz ~ 230   
% orid = 269; % seaz ~ 250   

% sta1 = 'J53C'; % close to trench
% sta1 = 'J43C'; % close to trench
% sta2 = 'J31C'; % close to ridge
% sta2 = 'J39C'; % close to ridge
% sta2 = 'J47C'; % close to ridge

%% ========== Load the two traces ==========
evt = dir([datdir,num2str(orid),'*']);
load([datdir,evt.name,'/_EQAR_',phase,'_',comp])
% load('synthetic_tstar/eg_eqar_456ST.mat')

ind1 = find(strcmp({eqar.sta},sta1)); %close to trench %  1  1  1  2
ind2 = find(strcmp({eqar.sta},sta2)); %close to ridge % 48 -3  5 30

dat1 = eqar(ind1).datT';
% tt10 = (eqar(ind1).tt'- eqar(ind1).pred_arrT)*spd;
tt1 = (eqar(ind1).tt'- eqar(ind1).abs_arrT)*spd + eqar(ind1).dT;
if ~isstr(ind2)
    dat2 = eqar(ind2).datT';
%     tt2 = (eqar(ind2).tt'- eqar(ind2).pred_arrT)*spd;
    tt2 = (eqar(ind2).tt'- eqar(ind2).abs_arrT)*spd + eqar(ind2).dT;
end

samprate = eqar(ind1).samprate;
dt = 1./samprate;
fnq = samprate/2;
T = prex+postx;
fmax = mean([eqar([ind1,ind2]).fcross]);

%interp to twin
tt = [-pretime:dt:pretime-dt]';
dat1 = interp1(tt1,dat1,tt);
dat2 = interp1(tt2,dat2,tt);

%% MAKE STACK INSTEAD
if ischar(ind2)
indgd = find(~isnan([eqar.dT]));
tt2 = [-pretime:dt:pretime]'; tt2 = tt2(1:length(tt1));
for ig = 1:length(indgd)
    is = indgd(ig);
    att = eqar(is).tt-eqar(is).abs_arrT; % shift to since MEASURED absolute arrival
    all_dat0(:,is) = interp1(att,eqar(is).datT,tt2,'linear',0);
end % loop on stas
indgd = 1:size(eqar);
indgd(mean(abs(all_dat0(:,indgd)))<1e-6)     = []; % kill zero traces
indgd(isnan(mean(abs(all_dat0(:,indgd))))) = []; % kill nan traces
stak = sum(all_dat0(:,indgd),2)/length(indgd);
dat2 = stak;
end

[ dat1 ] = filt_quick( dat1,1./40,10,dt);
[ dat2 ] = filt_quick( dat2,1./40,10,dt);

%% Plot two traces
figure(1), clf, set(gcf,'pos',[100 550 1000 350]), hold on
plot(tt,dat1./max(abs(dat1)),'k','LineWidth',1.5)
plot(tt,dat2./max(abs(dat2)),'r','LineWidth',1.5)
plot(-prex*[1 1],max(abs(dat1))*[-1 1],'b--',postx*[1 1],max(abs(dat1))*[-1 1],'b--')

%% Make set of period windows for bandpass filter
Tmids = logspace(log10(Tmin),log10(Tmax),Nwds)';
Twdhs = 0.5*diff(logspace(log10(Tmin/2),log10(2*Tmax),Nwds+1)');
% Twdhs = 0.5*linspace(1,3,Nwds);
fmids = 1./Tmids;



%% Test combspectra function
parms.comb.Tmin = Tmin; % fnq/2
parms.comb.Tmax = Tmax;
parms.comb.Nwds = Nwds;
parms.comb.Tw_opt = 'scale';
parms.comb.npol = 4;

parms.wind.pretime = pretime;
parms.wind.prex = prex;
parms.wind.postx = postx;
parms.wind.taperx = 0.1;

parms.qc.minacor = minacor;
parms.qc.maxphi = maxphi;
parms.qc.mtmRcomp = true; % option to compare comb and mtm estimates of As
parms.qc.focus_threshold = 3; % threshold value for focus metric to flag up as bad

parms.inv.amp2phiwt = amp2phiwt;
parms.inv.fmin = .1;
parms.inv.fmax = fmax;
parms.inv.ifwt = true;
parms.inv.corr_c_skip = true;
parms.inv.fmin_cskip = 0.1;
parms.inv.R2default = 0.9; % if zero, do not weight by R2 of each pairwise fit, Else, weight by multiple of this default.
parms.inv.opt = 1; %USE 1!   for calc_fdependent...  1 is all in one method; 2 is one-by-one method

parms.inv.alpha = alp;

% GET station-by-station acceptable f range from crossing freqs
fcross = [[eqar([ind1,ind2]).fcrosslo]',[eqar([ind1,ind2]).fcrosshi]'];
fcross(:,2) = 0.3;
ifplot = false;
% 
% tic
% [delta_tstar,delta_T,std_dtstar,pairwise] = combspectra([dat1,dat2],samprate,parms,1);
% fprintf('combspectra dtstar est\t= %.3f \n',diff(delta_tstar))
% fprintf('combspectra dT est\t= %.3f \n',diff(delta_T))
% toc
% tic
% [delta_tstar_skip,delta_T_skip,std_dtstar_skip,pairwise] = combspectra([dat1,dat2],fcross,samprate,eqar(1).par_dtstar_comb,1);
[delta_tstar_skip,delta_T_skip,std_dtstar_skip,pairwise,fmids] = combspectra([dat1,dat2],fcross,samprate,parms,combspectraplotopt);
dtstar_specR = eqar(ind2).dtstar - eqar(ind1).dtstar;
dtstar_comb  = eqar(ind2).dtstar_comb - eqar(ind1).dtstar_comb;

fprintf('whole event dtstar_specR est\t= %.3f \n',dtstar_specR)
fprintf('whole event dtstar_comb est\t= %.3f \n',dtstar_comb)
fprintf('whole event dT est    \t= %.3f \n',eqar(ind2).dT-eqar(ind1).dT)
fprintf('BETTER? combspectra dtstar est\t= %.3f \n',diff(delta_tstar_skip))
fprintf('BETTER? combspectra dT est\t= %.3f \n',diff(delta_T_skip))
% toc

[gdofgd,focusscore] = amp2phi_QC_focusing(parms,pairwise,fmids,1,{sta1,sta2});
fprintf('Focus score:\n')
disp(focusscore)

%% Get ready to plot
As = pairwise.As';
phis = pairwise.phis';
wts = pairwise.wts';
inds = find(fmids<=fmax & abs(phis)<maxphi & sqrt(wts)>minacor);
[ dtstar,dT,A0,misfit,res ] = invert_1pair_Aphi_4_dtdtstar( As(inds),phis(inds),fmids(inds), wts(inds),amp2phiwt,alp);
frq = eqar(ind1).frq;

%% NEW - diff. phase and amplitude spectra from fourier decomp. 
% subset data to signal windows
cp = struct('samprate',1./dt,'pretime',pretime,'prex',prex,'postx',postx,...
                'taperx',0.1,'fhi',100,'flo',0,'npoles',2,'norm',0);
[ dat1w ] = data_clean( dat1,cp); 
[ dat2w ] = data_clean( dat2,cp);

[dat1_fft,ff_fft] = fft_ze(dat1w,dt);
[dat2_fft,~]   = fft_ze(dat2w,dt);
N = length(dat1_fft);
fplt = (ff_fft>=max(fcross(:,1))) & (ff_fft<=min(fcross(:,2)));
specRfft = smooth(log(abs(dat2_fft)./abs(dat1_fft)));
% specPfft = smooth(unwrap(angle(dat2_fft./dat1_fft))./(-2*pi*ff_fft));
specPfft = smooth((angle(dat2_fft./dat1_fft))./(-2*pi*ff_fft));

specRfft = specRfft(fplt);
specPfft = specPfft(fplt);
ff_fft = ff_fft(fplt);
% wt_fft = 1./(0.5+abs(interp1(eqar(ind1).frq,log(eqar(ind2).specs./eqar(ind1).specs),ff_fft) - specRfft));
wt_fft_a = 0.4./(.2+abs(interp1(eqar(ind1).frq,log(eqar(ind2).specs./eqar(ind1).specs),ff_fft) - specRfft));
wt_fft_p = [1;abs(diff(specPfft)./diff(ff_fft))<10];
wt_fft = (wt_fft_a.*wt_fft_p) + 0.0001;
figure(444),clf
subplot(211), hold on
plot(ff_fft,specRfft ,'-b','linewidth',2)
plot(eqar(ind1).frq,log(eqar(ind2).specs./eqar(ind1).specs),'-r')
plot(ff_fft,wt_fft,'k')
xlim([0 0.5])
subplot(212), hold on
plot(ff_fft,specPfft ,'-b','linewidth',2)
xlim([0 0.5])




%% best-fit f dependence
alpha_test = [0:0.05:0.9]';
for ia = 1:length(alpha_test)
    [ alpha_dtstar(ia),~,~,alpha_misfit(ia),~ ] = invert_1pair_Aphi_4_dtdtstar( As(inds),phis(inds),fmids(inds), wts(inds),amp2phiwt,alpha_test(ia));
end
figure(11), clf
plot(alpha_test,alpha_misfit,'o-')
figure(12), clf
plot(alpha_test,alpha_dtstar,'o-')

%% predicted slopes from this pair of stas
ws = 2*pi*fmids;
if alp == 0
    lnApred   = log(A0) - pi*fmids*dtstar;
    phipred = -(1/pi)*log(fmids)*dtstar + dtstar + dT;
else
    lnApred = log(A0) - 0.5 * ws.^(1-alp) * w0^alp * dtstar;
    phipred = 0.5*ws.^(-alp).*(2*pi).^alp * cot(alp*pi/2)*dtstar + dT;
end

%% predicted slopes from inversion of all stas - in eqar
% specR
lnApred_specR   = log(A0) - pi*fmids*dtstar_specR;
phipred_specR = -(1/pi)*log(fmids)*dtstar_specR + dtstar_specR + dT;

% comb
lnApred_comb   = log(A0) - pi*fmids*dtstar_comb;
phipred_comb = -(1/pi)*log(fmids)*dtstar_comb + dtstar_comb + dT;


%% do for alp0.27 too
[ dtstar_a27,dT_a27,A0_a27,misfit,res ] = invert_1pair_Aphi_4_dtdtstar( As(inds),phis(inds),fmids(inds), wts(inds),amp2phiwt,0.27);
lnApred_a27 = log(A0_a27) - 0.5 * ws.^(1-0.27) * w0^0.27 * dtstar_a27;
phipred_a27 = 0.5*ws.^(-0.27).*(2*pi).^0.27 * cot(0.27*pi/2)*dtstar_a27 + dT_a27;

%% do for fft spectra too
[ dtstar_fft,dT_fft,A0_fft] = invert_1pair_Aphi_4_dtdtstar( exp(specRfft),specPfft,ff_fft, wt_fft,amp2phiwt,0);
lnApred_fft = log(A0_fft) - pi*ff_fft * dtstar_fft;
phipred_fft = -(1/pi)*log(ff_fft)*dtstar_fft + dtstar_fft + dT_fft;
fprintf('BETTER? specXPR dtstar est = %.3f \n',dtstar_fft)
fprintf('BETTER? specXPR dT est\t   = %.3f \n',dT_fft)


    
%% =========================================================================
%% ====================== START THE PLOTS ==================================

figure(30), clf, set(gcf,'pos',[600 309 736 1036])

% ================ PLOT THE WAVEFORMS ================
subplot(5,1,1), hold on
plot(tt,dat1./max(abs(dat1)),'k','LineWidth',2)
plot(tt,dat2./max(abs(dat2)),'r','LineWidth',2)
plot(-prex*[1 1],[-1 1],'b--',postx*[1 1],[-1 1],'b--','linewidth',1.5)
text(-prex-18,-0.5,sprintf('\\textbf{%s}',sta1),'interpreter','latex','fontsize',12,'horizontalalignment','left','verticalalignment','bottom','color','k')
text(-prex-18,-0.8,sprintf('\\textbf{%s}',sta2),'interpreter','latex','fontsize',12,'horizontalalignment','left','verticalalignment','bottom','color','r')
text(postx+15,-1.65,'Time $\rightarrow$','interpreter','latex','fontsize',15,'horizontalalignment','left','verticalalignment','bottom')


% xlabel('Time, s','FontSize',22,'interpreter','latex')
set(gca,'xlim',[-prex-20 postx+20],'ylim',[-1.1 1.1],...%'XAxisLocation','top',...
    'Fontsize',12,'linewidth',2,'box','on','ytick',[])

% ================ PLOT THE AMPLITUDE SPECTRA ================
subplot(5,1,2:3), hold on
% plot(frq(2:end),log(spec2./spec1),'-ok')
plot(eqar(ind1).frq,log(eqar(ind2).specs./eqar(ind1).specs),'-xk','linewidth',1.5)
% plot(eqar(ind1).frq,log(eqar(ind2).specs./eqar(ind2).specn),'-ob')
% plot(eqar(ind1).frq,log(eqar(ind1).specs./eqar(ind1).specn),'-om')
plot(fmids,lnApred,'g','Linewidth',2.5)
plot(fmids,lnApred_a27,'--g','Linewidth',1.5)
% plot predictions from all evts...
plot(fmids,lnApred_specR,'r','Linewidth',1.5)
plot(fmids,lnApred_comb,'c','Linewidth',1.5)

plot(ff_fft,specRfft ,'-b','linewidth',2)
scatter(ff_fft,specRfft,40*wt_fft,'markerfacecolor','b')
plot(ff_fft,lnApred_fft,'--b','Linewidth',1.5)

gdf = wts~=0;
scatter(fmids(gdf),log(As(gdf)),150*wts(gdf),'or','MarkerFaceColor','r')
scatter(fmids(inds),log(As(inds)),150*wts(inds),'ok','MarkerFaceColor','r','linewidth',1.5)

% ha1 = plot([eqar(ind1).fcross,1],[eqar(ind1).fcross,0],1.9,90,4,0.12,'FaceColor','k'); % fcross for sta 1
% ha2 = arrow([eqar(ind2).fcross,1],[eqar(ind2).fcross,0],1.9,90,4,0.12,'FaceColor','k'); % fcross for sta 2
% text(eqar(ind1).fcross,1.3,'$f^{max}$','interpreter','latex','fontsize',18,'horizontalalignment','center','verticalalignment','bottom')
% text(eqar(ind2).fcross,1.3,'$f^{max}_2$','interpreter','latex','fontsize',16,'horizontalalignment','center','verticalalignment','bottom')
plot(max(fcross(:,1))*[1 1],[-1 1],'--b','LineWidth',1.5)
plot(min(fcross(:,2))*[1 1],[-1 1],'--b','LineWidth',1.5)

% xlabel('freq (Hz)','FontSize',16,'interpreter','latex')
ylabel('$\ln\,(R_{12})$','FontSize',22,'interpreter','latex')
set(gca,'fontsize',15,'Xscale','linear','xticklabel',[],'xlim',[0.00 0.6],'ylim',[-3.5 2.5],'linewidth',2,'box','on')

% ================ PLOT THE PHASE SPECTRA ================
subplot(5,1,4:5), hold on
plot(fmids,phipred,'g','Linewidth',2.5)
plot(fmids,phipred_a27,'--g','Linewidth',1.5)
% plot predictions from all evts...
plot(fmids,phipred_specR,'r','Linewidth',1.5)
plot(fmids,phipred_comb,'c','Linewidth',1.5)

plot(ff_fft,specPfft ,'-b','linewidth',2)
scatter(ff_fft,specPfft,40*wt_fft,'markerfacecolor','b')
plot(ff_fft,phipred_fft ,'--b','linewidth',2)

scatter(fmids(gdf),phis(gdf),150*wts(gdf),'or','MarkerFaceColor','r')
scatter(fmids(inds),phis(inds),150*wts(inds),'ok','MarkerFaceColor','r','linewidth',1.5)

% ha1 = arrow([eqar(ind1).fcross,1],[eqar(ind1).fcross,0],5,90,2,0.1,'FaceColor','m'); % fcross for sta 1
% ha2 = arrow([eqar(ind2).fcross,1],[eqar(ind2).fcross,0],5,90,2,0.1,'FaceColor','g'); % fcross for sta 2
% plot(fmax*[1 1],[-1 1],'--b')  
text(0.5,-.2,['$\mathbf{\Delta t^* = ',num2str(dtstar,'%.1f'),'}$'],'interpreter','latex','fontsize',20,'horizontalalignment','right','verticalalignment','bottom')

text(0.03,0,'More delayed $\rightarrow$','interpreter','latex','fontsize',16,...
    'rotation',90,'horizontalalignment','left','verticalalignment','bottom')


xlabel('frequency (Hz)','FontSize',22,'interpreter','latex')
ylabel('$\Delta \psi_{12}$ (s)','FontSize',22,'interpreter','latex')
set(gca,'fontsize',15,'Xscale','linear','xlim',[0.00 0.6],'linewidth',2,'box','on')


% ================ SAVE FIGURE ================
if ifsave
    ofile = sprintf('eg_comb_orid%.0f_%s%s_%s_v_%s_alp%.2f',orid,phase,comp,sta1,sta2,alp);
    save2pdf(30,ofile,'figs');
end