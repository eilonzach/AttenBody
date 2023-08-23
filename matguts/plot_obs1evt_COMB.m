function eqar = plot_obs1evt_COMB( eqar,pairwise)

datwind = [-70 100];
filtfs = [0.05 2]; % in hz
tstlim = 1.5*[-1 1];
agelim = [0 10];
% refsta = 'LON';
amp_seis = 1.5; % amplification factor for the traces

t = tauptime('p',eqar(1).phase,'deg',eqar(1).gcarc); 
evtime = eqar(1).pred_arrT - t(1).time;

% only use "good" stations
gd=[]; 
for is = 1:length(eqar), 
    if ~isempty(eqar(is).dtstar_comb) && ~isnan(eqar(is).dtstar_comb), 
        gd = [gd;is]; 
    end
end
eqar = eqar(gd);
% put NaNs in for stations without specR dtstar:
for is=1:length(eqar), if isempty(eqar(is).dtstar), eqar(is).dtstar=NaN; end; end

%% get deets
phase = eqar(1).phase;
comp = eqar(1).par_dtstar_specR.comp;
window = [-eqar(1).par_dtstar_comb.wind.prex,eqar(1).par_dtstar_comb.wind.postx];
fmax = nanmean([eqar.fcross]');
if isempty(filtfs), filtfs = eqar(1).par_dtstar_specR.filtfs; end

%% Find the average spectrum and use that as a reference
specss_ref = zeros(size(eqar(1).specss));
kk = 0;
for is = 1:length(eqar)
    if isempty(eqar(is).specss), continue; end
    specss_ref = specss_ref + eqar(is).specss;
    kk = kk+1;
end
specss_ref = specss_ref/kk;

%% ===================================================================
%% ------------------ MAP WITH DTSTAR FOR THIS EVENT -----------------
%% ===================================================================

figure(32), clf, hold on
mkfig_CascMAP
set(gca,'xlim',[-134.5,-120],'ylim',[42,50],'fontsize',24,...
    'xtick',[-135:5:-120],'ytick',[40:2:50])
scatter([eqar.slon],[eqar.slat],300,[eqar.dtstar_comb],'filled','MarkerEdgeColor','k')
% % see where the slab depth section is
% scatter(slabdepth(:,3),slabdepth(:,4),10,.05*slabdepth(:,2))

% text([eqar.slon]+0.2,[eqar.slat],{eqar.sta})
cmap = parula;
colormap(cmap)
caxis(tstlim)

% plot section
plot(section_lola(:,1),section_lola(:,2),'k','LineWidth',2)
plot(linterp(section_x,section_lola(:,1),[-100:100:700]'),...
     linterp(section_x,section_lola(:,2),[-100:100:700]'),...
     '.k','MarkerSize',25)
text(linterp(section_x,section_lola(:,1),[-100:200:700]'),...
     linterp(section_x,section_lola(:,2),[-100:200:700]')+0.14,...
     num2str([-100:200:700]'),'FontWeight','bold','fontsize',15)

%% colour bar
cbar_custom(gca, 'location',[-133.5 -133.1 42.5 45.5],'tickside','right',...
    'lims',tstlim,'tickvals',[tstlim(1):0.5:tstlim(2)],'cmap',cmap,...
    'FontSize',17,'FontWeight','bold',...
	'title',sprintf('$\\Delta t^*_%s$ \\,(s)',eqar(1).phase),'interpreter','latex');

%% ===================================================================
%% -------------------------- DO SOME CALCS --------------------------
%% ===================================================================

%% work out distance to ridge
Xrdg = dist_to_rdg([eqar.slat]',[eqar.slon]'); 
eqar(1).Xrdg = Xrdg(1);
eqar = dealto(eqar,'Xrdg',Xrdg);
[~,Xord] = sort([eqar.Xrdg]); % sort by longitude
[~,Xord] = sort(Xord);
Xrdglims = [min(abs([eqar.Xrdg])),max([eqar.Xrdg])];
kpd = 78; % km per degree of longitude

%% load sed thickness
sectseds = dlmread('mapdata/2D_bathym_section_sediments.txt','\t',2,0); % load sediment thickness for section
sectseds = sectseds(:,[1,2,4]);

%% load + calc. slab depths
slabdepth = dlmread('/Users/zeilon/Work/CASCADIA/plots/depth_to_slab/slabgeom.txt',' ');
% prof_centre=-124.18427/47.09777, az = 100
[slabdepth(:,4),slabdepth(:,3)] = reckon_km(47.09777,-124.18427,slabdepth(:,1),100);
slabdepth(:,1) = slabdepth(:,1)-slabdepth(1,1)+250;
load('mapdata/2D_section_fakeslabdepth.mat'); slbz = faked_slab_depth;

%% Calc differential travel times just from bathymetry + seds
z_seflr = section_z;
sedZ = sectseds(:,3);
sedZ(section_x>320) = 1e3;
z_sedbo = z_seflr - sedZ;
% Make sectseds at least 1km thick on the shelf
% look at flueh
z_slbto = linterp(faked_slab_depth(:,1),faked_slab_depth(:,2),section_x);
z_ocmoh = z_slbto - 6.5e3;
z_fifty = -50e3*ones(size(z_seflr));

% sta_z_seflr = linterp(section_x,z_seflr,Xrdg');
% sta_z_sedbo = linterp(section_x,z_sedbo,Xrdg');
% sta_z_slbto = linterp(section_x,z_slbto,Xrdg');
% sta_z_ocmoh = linterp(section_x,z_ocmoh,Xrdg');
% sta_z_fifty = -50e3*ones(size(z_seflr));

v_asth = 4.4e3;
v_oc = 3.9e3;
v_cr = 3.6e3;
v_sed = 1.5e3;

TT = (z_seflr-z_sedbo)/v_sed + ...
     (z_sedbo-z_slbto)/v_cr + ...
     (z_slbto-z_ocmoh)/v_oc + ...
     (z_ocmoh-z_fifty)/v_asth;
TT = TT - nanmean(TT);

%% ===================================================================
%% --------------- SECTION WITH DTSTAR FOR THIS EVENT ----------------
%% ===================================================================
figure(17), clf, set(gcf,'pos',[59   258   800   913])
%% topo
maxy = 3000;
miny = -7000;
subplot(3,1,1), hold on
% fill in some layers
fill([section_x;flipud(section_x)],[zeros(size(section_x));flipud(z_seflr)],[221 253 251]/255,'EdgeColor','none')
fill([section_x;flipud(section_x)],[z_seflr;flipud(z_sedbo)],[253 240 221]/255,'EdgeColor','none')
fill([section_x;flipud(section_x)],[z_sedbo;flipud(z_slbto)],[255 241 236]/255,'EdgeColor','none')
fill([section_x;flipud(section_x)],[z_slbto;-1e5*ones(size(section_x))],[245 255 236]/255,'EdgeColor','none')
% draw bathym, seds, slab
plot(slbz(:,1),slbz(:,2),'--k','LineWidth',1)
plot(section_x,z_sedbo,':k','LineWidth',1)
plot(section_x,z_seflr,'k','LineWidth',2)


% figure things
xlim([-80 460]),ylim([miny maxy]+[-1 1])
set(gca,'visible','off','fontsize',12)
% draw axes
plot([-80 720],[0 0],'--k') % x-axis
plot([-80 -80],[miny maxy],'k','LineWidth',2)
plot([-80 -72],[miny miny],'k',[-80 -72],[maxy maxy],'k','LineWidth',1)
text(-80,maxy+1e3,num2str(maxy),'Fontsize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
text(-80,miny-1e3,num2str(miny),'Fontsize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
text(-115,mean([miny maxy]),'\textbf{Elev. (m)}','Fontsize',17,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','Rotation',90)
% title
text(180,miny-2e3,sprintf('%s $~$  $%s$-wave, %s-comp',epoch2str(evtime,'%Y-%m-%d %H:%M:%S'),phase,comp),'Fontsize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex')
% ridge axis label
plot(0,2500,'vk','MarkerSize',8,'MarkerFaceColor','k')
text(0,3500,'Axis','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',16,'Interpreter','latex')
% deformation front label
plot(305,2500,'vk','MarkerSize',8,'MarkerFaceColor','k')
text(305,3500,'DF','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',16,'Interpreter','latex')
% Coastline label
plot(440,2500,'vk','MarkerSize',8,'MarkerFaceColor','k')
text(440,3500,'Coastline','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',16,'Interpreter','latex')
% Arc label
plot(616,2500,'vk','MarkerSize',8,'MarkerFaceColor','k')
text(616,3500,'Arc','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',16,'Interpreter','latex')
% Slab top label
text(327,-5800,'Top of slab','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',15,'Interpreter','latex','rotation',-24)
% Stations
plot(Xrdg*kpd,[eqar.selev],'^b','MarkerSize',8,'MarkerFaceColor','r')

% delta tstar
subplot(3,1,2), hold on
scatter(Xrdg*kpd,[eqar.dtstar],20,'g','MarkerFaceColor','g')
scatter(Xrdg*kpd,[eqar.dtstar_comb],30./[eqar.stds_comb],'b','MarkerFaceColor','b')
% plot(Xrdg(ilan)*kpd,[eqar(ilan).dtstar],'.r','MarkerSize',18)
% figure things
grid off
set(gca,'xlim',[-80 460],'ylim',[-3 3],'Fontsize',17,'xticklabel',[])
ax = axis;
plot(ax(1:2),[0 0],':k','Linewidth',1.5)
text(-115,0,'$\mathbf{\Delta t^{\ast}_S}$','Fontsize',20,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','Rotation',90)

% delta t
subplot(3,1,3), hold on
plot(section_x,TT + 1.20,'--r','LineWidth',1.5) % < plot on the prediction of TT from just topog + seds.
scatter(Xrdg*kpd,[eqar.dT],80,'b','MarkerFaceColor','b')
% plot(Xrdg*kpd,[eqar.dT_comb],'.b','MarkerSize',18)
% plot(Xrdg(ilan)*kpd,[eqar(ilan).dT],'.r','MarkerSize',18)
% figure things
grid off
set(gca,'xlim',[-80 460],'ylim',[-2.5 2.5],'Fontsize',17)
ax = axis;
plot(ax(1:2),[0 0],':k','Linewidth',1.5)
text(-115,0,'$\mathbf{\delta T_S}$','Fontsize',20,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','Rotation',90)
text(mean(ax(1:2)),-3.5,'\textbf{Distance from ridge (km)}','Fontsize',18,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex')

%% ===================================================================
%% ---------------------------- WAVEFORMS ----------------------------
%% ===================================================================

for is = 1:length(eqar)
    dat = eqar(is).(['dat',comp])';

    W = 2*filtfs./eqar(is).samprate;
    [fb, fa] = butter(2,W);
    dd(:,is) = filtfilt(fb,fa,dat);
end

%    plot
figure(2), clf, set(gcf,'position',[150 000 600 800]), hold on
for is = 1:length(eqar)
	tt = eqar(is).tt-eqar(is).abs_arrT;

    hp = plot(tt,amp_seis*dd(:,is)/max(max(abs(dd))) + Xord(is),'LineWidth',2);
    set(hp,'color',colour_get(eqar(is).staage,agelim(2),agelim(1),flipud(jet)))
    set(hp,'color',colour_get(abs(eqar(is).Xrdg),400/kpd,0,flipud(jet)))
%     text(datwind(1)+2,Xord(is)+0.4,eqar(is).sta,...                        % station names
%         'FontSize',8,'interpreter','latex','HorizontalAlignment','left')
%     text(datwind(2) + 1,Xord(is),sprintf('%.0f km',eqar(is).Xrdg*kpd),...  % distance to ridge
%         'FontSize',10,'interpreter','latex')
end
plot([1;1]*window,[0;length(eqar)+1]*[1 1],'--k','LineWidth',1)
set(gca,'FontSize',15,'YTick',[],'position',[0.16 0.11 0.75 0.815])
axis([datwind,-0.5,length(eqar)+1.5])
xlabel('\textbf{Time from phase onset (s)}','FontSize',18,'interpreter','latex')
title(sprintf('%s $~$ %s-wave, %s-comp, $f_{hi}~$: %.2f',epoch2str(evtime,'%Y-%m-%d %H:%M:%S'),...
    phase,comp,fmax),'Fontsize',16,'Interpreter','latex')

% ridge axis label
text(datwind(1)-2,sum([eqar.Xrdg] < 0)+0.5,'\textbf{Axis $\succ$}',...
    'HorizontalAlignment','right','VerticalAlignment','bottom',...
    'FontSize',16,'FontWeight','bold','Interpreter','latex')
% deformation front label
text(datwind(1)-2,sum([eqar.Xrdg] < 305/kpd)+0.5,'\textbf{DF $\succ$}',...
    'HorizontalAlignment','right','VerticalAlignment','bottom',...
    'FontSize',16,'FontWeight','bold','Interpreter','latex')
% Coastline label
text(datwind(1)-2,sum([eqar.Xrdg] < 440/kpd)+0.5,'\textbf{Coast $\succ$}',...
    'HorizontalAlignment','right','VerticalAlignment','bottom',...
    'FontSize',16,'FontWeight','bold','Interpreter','latex')
% Arc label
% text(-32,sum([eqar_obs.Xrdg] < 616/kpd)+0.5,'Arc $\rightarrow$',...
%     'HorizontalAlignment','right','VerticalAlignment','bottom',...
%     'FontSize',15,'FontWeight','bold','Interpreter','latex')


%% ===================================================================
%% --------------------------- SPECTRA  ------------------------------
%% ===================================================================
if nargin > 1 && ~isempty(pairwise)
%% spectra
figure(88),clf, set(gcf,'position',[400 100 580 700])

age_bins = [0:1:8];
fmids = pairwise.fmids;
Nb = length(age_bins);
Nf = length(fmids);
As_bin = zeros(Nb,Nf);
phis_bin = zeros(Nb,Nf);
wts_bin = zeros(Nb,Nf);

for ib = 1:Nb-1
    indbin = ([eqar.staage] <= age_bins(ib+1)) & ([eqar.staage] >= age_bins(ib))
    As_bin(ib,:) = nanmean(pairwise.As(find(indbin),:),1); % do actually want 'find', because only first Nstas of pairwise are relevant for the stack
    phis_bin(ib,:) = nanmean(pairwise.phis(find(indbin),:),1);
    wts_bin(ib,:) = nanmean(pairwise.wts(find(indbin),:),1);
end    

ax1 = axes('position',[0.13  0.54  0.775  0.38]); hold on
ax2 = axes('position',[0.13  0.11  0.775  0.38]); hold on

for ib = 1:Nb-1
    scatter(ax1,fmids,log(As_bin(ib,:)),130*wts_bin(ib,:),'ok','Linewidth',0.5,'MarkerFaceColor',...
            colour_get(mean(age_bins(ib:ib+1)),agelim(2),agelim(1),flipud(jet)))
    plot(ax1,fmids,log(As_bin(ib,:)),'--','Linewidth',2,'color',...
            colour_get(mean(age_bins(ib:ib+1)),agelim(2),agelim(1),flipud(jet)))
    set(ax1,'Xscale','linear','xlim',[0.03 fmax+0.05],'ylim',[-1.5 2.5],'xticklabel',[],'ytick',[-1.5:1:2.5],...
        'fontsize',15,'box','on','linewidth',2)
    ylabel(ax1,'$\mathbf{\ln {(R)}}$','FontSize',18,'Interpreter','latex')

    scatter(ax2,fmids,phis_bin(ib,:),130*wts_bin(ib,:),'ok','Linewidth',0.5,'MarkerFaceColor',...
            colour_get(mean(age_bins(ib:ib+1)),agelim(2),agelim(1),flipud(jet)))
    plot(ax2,fmids,phis_bin(ib,:),'--','Linewidth',2,'color',...
            colour_get(mean(age_bins(ib:ib+1)),agelim(2),agelim(1),flipud(jet)))
    set(ax2,'Xscale','linear','xlim',[0.03 fmax+0.05],'ylim',[-3 3],'fontsize',15,'box','on','linewidth',2)
    xlabel(ax2,'\textbf{frequency (Hz)}','FontSize',18,'Interpreter','latex'), 
    ylabel(ax2,'$\mathbf{\Delta \psi}$ \,\,\,\textbf{(s)}','FontSize',18,'Interpreter','latex')
end

% cbar_custom(gca, 'location',[0.06 .12 1 1.18],'tickside','bottom',...
%     'lims',kpd*Xrdglims_plot,'tickvals',round_level([kpd*Xrdglims_plot(1):50:kpd*Xrdglims_plot(2)],100),...
%     'FontSize',10,'FontWeight','bold','cmap',flipud(jet),...
%     'title','Distance E of ridge (km)','interpreter','latex');
% 



% for is = 1:length(eqar)
%     subplot(211), hold on
%     scatter(fmids,log(pairwise.As(is,:)),110*pairwise.wts(is,:),'ok','Linewidth',0.5,'MarkerFaceColor',...
%             colour_get(eqar(is).staage,12,0,flipud(jet)))
%     plot(fmids,log(pairwise.As(is,:)),...
%             'color',colour_get(eqar(is).staage,max([eqar.staage]),min([eqar.staage]),flipud(jet)))
%     set(gca,'Xscale','linear','xlim',[0.03 fmax+0.1],'ylim',[-3.5 3.5],'fontsize',14,'box','on','linewidth',2)
%     %     xlabel('\textbf{frequency (Hz)}','FontSize',18,'Interpreter','latex')
%     ylabel('$\mathbf{\ln {(R)}}$','FontSize',18,'Interpreter','latex')
% 
%     subplot(212), hold on
%     scatter(fmids,pairwise.phis(is,:),110*pairwise.wts(is,:),'ok','Linewidth',0.5,'MarkerFaceColor',...
%             colour_get(eqar(is).staage,max([eqar.staage]),min([eqar.staage]),flipud(jet)))
%     plot(fmids,pairwise.phis(is,:),':','color',...
%             colour_get(eqar(is).staage,12,0,flipud(jet)))
%     set(gca,'Xscale','linear','xlim',[0.03 fmax+0.1],'ylim',[-3.5 4.5],'fontsize',14,'box','on','linewidth',2)
%     xlabel('\textbf{frequency (Hz)}','FontSize',18,'Interpreter','latex'), 
%     ylabel('$\mathbf{\Delta \psi}$ \,\,\,\textbf{(s)}','FontSize',18,'Interpreter','latex')
% end



end

end

