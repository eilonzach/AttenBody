function plot_ATTEN_TandF_domain_COMB( eqar, plotopt )
% plotopt can be "all" or "obs"

if nargin<2
    plotopt = 'all';
end

datwind = [-70 100];
filtfs = [0.05 2]; % in hz
tstlim = 1.5*[-1 1];
% refsta = 'LON';

kpd = 110; % kilometers per degree
spd = 24*60*60; % sec per day, to convert from serial time

t = tauptime('p',eqar(1).phase,'deg',eqar(1).gcarc); 
evtime = eqar(1).pred_arrT - t(1).time;

% get "good" stations
gd=[]; for is = 1:length(eqar), if ~isempty(eqar(is).dtstar_comb) && ~isnan(eqar(is).dtstar_comb), gd = [gd;is]; end, end
eqar_gd = eqar(gd);
% put NaNs in for stations without specR dtstar:
for is=1:length(eqar_gd), if isempty(eqar_gd(is).dtstar), eqar_gd(is).dtstar=NaN; end; end

%% get deets
phase = eqar_gd(1).phase;
comp = eqar_gd(1).par_dtstar_specR.comp;
window = [-eqar_gd(1).par_dtstar_comb.wind.prex,eqar_gd(1).par_dtstar_comb.wind.postx];
fmax = eqar_gd(1).par_dtstar_comb.inv.fmax;
if isempty(filtfs), filtfs = eqar_gd(1).par_dtstar_specR.filtfs; end

% % << MAKE DEFUNCT>>
% iobs = ~cellfun('isempty',regexp({eqar_gd.sta},'([J,F,M,G]*[0-9][0-9][A-C]$)'));
% ilan =  cellfun('isempty',regexp({eqar_gd.sta},'([J,F,M,G]*[0-9][0-9][A-C]$)'));
% 
% % account for stupidly named land stations
% ilan(find(strcmp({eqar_gd(iobs).sta},'M02C') & [eqar_gd(iobs).slat]==41.392))=true;
% iobs(find(strcmp({eqar_gd(iobs).sta},'M02C') & [eqar_gd(iobs).slat]==41.392))=false;


% if nargin<2 || isempty(refsta)
% mlat = mean([eqar_gd.slat]); 
% mlon = mean([eqar_gd.slon]); 
% dlalo = distance(mlat,mlon,[eqar_gd.slat],[eqar_gd.slon]);
% iref = find(dlalo==min(dlalo));
% else
% iref = find(strcmp({eqar_gd.sta},refsta));
% end

%% Find the average spectrum and use that as a reference
specss_ref = zeros(size(eqar_gd(1).specss));
kk = 0;
for is = 1:length(eqar_gd)
    if isempty(eqar_gd(is).specss), continue; end
    specss_ref = specss_ref + eqar_gd(is).specss;
    kk = kk+1;
end
specss_ref = specss_ref/kk;


%% Sort by distance from MER
% plotting things
addpath('PLOTTING')
map_parameters

Xrft = dist2line([MER.lon(1),MER.lat(1)],[MER.lon(2),MER.lat(2)],[[eqar_gd.slon]' [eqar_gd.slat]']);
Xrft = distance(34.5,11,[eqar_gd.slon]',[eqar_gd.slat]');
eqar_gd = dealto(eqar_gd,'Xrft',Xrft);
[~,Xord] = sort([eqar_gd.Xrft]); % sort by longitude
[~,Xord] = sort(Xord);
Xrftlims_plot = [0 ceil(max(abs(Xrft)))]; if Xrftlims_plot(2)==0; Xrftlims_plot(2) = 1e-5; end

%% ===================================================================
%% ----------------- AGREEMENT BETWEEN COMB AND SPECR ----------------
%% ===================================================================
figure(44);hold on; 
for is = 1:length(eqar_gd)
    plot(eqar(is).dtstar_comb,eqar(is).dtstar,'ok');
end
xlabel('Comb'); ylabel('specR');
 
%% ===================================================================
%% ------------------ MAP WITH DTSTAR FOR THIS EVENT -----------------
%% ===================================================================

figure(32), clf, hold on
% mkfig_EARmap
axBM = gca;

% new axes
axAT = axes('pos',get(axBM,'pos'));
m_coord('geographic');

m_scatter([eqar_gd.slon]',[eqar_gd.slat]',200,[eqar_gd.dtstar_comb]','filled')
m_text([eqar_gd.slon]'+0.1,[eqar_gd.slat]'-0.1,{eqar_gd.sta},'fontsize',10,'fontweight','bold')
m_grid('backcolor','none','fontsize',17,'box','on','linewidth',2.5); pause(0.01);
% text([eqar_gd.slon]+0.2,[eqar_gd.slat],{eqar_gd.sta})

cmap = parula;
colormap(axAT,cmap)
caxis(axAT,tstlim)

% % plot section
% plot(section_lola(:,1),section_lola(:,2),'k','LineWidth',2)
% plot(linterp(section_x,section_lola(:,1),[-100:100:700]'),...
%      linterp(section_x,section_lola(:,2),[-100:100:700]'),...
%      '.k','MarkerSize',25)
% text(linterp(section_x,section_lola(:,1),[-100:100:700]'),...
%      linterp(section_x,section_lola(:,2),[-100:100:700]')+0.12,...
%      num2str([-100:100:700]'),'FontWeight','bold')

%% colour bar
[cblo,cbla]=m_ll2xy([33;36.5],[15;15.5]);
cbar_custom(axAT, 'location',[cblo' cbla'],'tickside','bottom',...
    'lims',tstlim,'tickvals',[tstlim(1):0.5:tstlim(2)],'cmap',cmap,...
    'FontSize',12,'FontWeight','bold',...
	'title',sprintf('$\\Delta t^*_%s$ \\,(s)',eqar(1).phase),'interpreter','latex');
 


% %% ===================================================================
% %% --------------- SECTION WITH DTSTAR FOR THIS EVENT ----------------
% %% ===================================================================
% figure(17), clf, set(gcf,'pos',[59   258   871   613])
% %% topo
% subplot(5,1,1), hold on
% plot(section_x,section_z,'k','LineWidth',2)
% % figure things
% xlim([-80 720]),ylim([-5001 5001])
% set(gca,'visible','off','fontsize',12)
% % draw axes
% plot([-80 720],[0 0],'--k') % x-axis
% plot([-80 -80],[-5000 5000],'k','LineWidth',2)
% plot([-80 -72],[-5000 -5000],'k',[-80 -72],[5000 5000],'k','LineWidth',1)
% text(-80,7000,'5000','Fontsize',10,'HorizontalAlignment','center','VerticalAlignment','middle')
% text(-80,-7000,'-5000','Fontsize',10,'HorizontalAlignment','center','VerticalAlignment','middle')
% text(-115,0,'Elev (m)','Fontsize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','Rotation',90)
% % title
% text(300,-7000,sprintf('%s $~$  $%s$-wave, %s-comp',epoch2str(evtime,'%Y-%m-%d %H:%M:%S'),phase,comp),'Fontsize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex')
% % ridge axis label
% plot(0,4500,'vk','MarkerSize',8,'MarkerFaceColor','k')
% text(0,6500,'Axis','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'Interpreter','latex')
% % deformation front label
% plot(305,4500,'vk','MarkerSize',8,'MarkerFaceColor','k')
% text(305,6500,'DF','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'Interpreter','latex')
% % Coastline label
% plot(440,4500,'vk','MarkerSize',8,'MarkerFaceColor','k')
% text(440,6500,'Coastline','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'Interpreter','latex')
% % Arc label
% plot(616,4500,'vk','MarkerSize',8,'MarkerFaceColor','k')
% text(616,6500,'Arc','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'Interpreter','latex')
% 
% 
% % delta tstar
% subplot(5,1,2:3), set(gca,'fontsize',12), hold on
% plot(Xord(iobs)*kpd,[eqar_gd(iobs).dtstar],'.g','MarkerSize',10)
% plot(Xord(iobs)*kpd,[eqar_gd(iobs).dtstar_comb],'.b','MarkerSize',18)
% plot(([eqar_gd(ilan).slon]+124)*kpd + 440,[eqar_gd(ilan).dtstar],'.g','MarkerSize',10)
% plot(([eqar_gd(ilan).slon]+124)*kpd + 440,[eqar_gd(ilan).dtstar_comb],'.r','MarkerSize',18)
% % plot(Xord(ilan)*kpd,[eqar_gd(ilan).dtstar],'.r','MarkerSize',18)
% % figure things
% grid on
% text(-115,0.5,'$\Delta t^{\ast}_S$','Fontsize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','Rotation',90)
% xlim([-80 720])
% 
% % delta t
% subplot(5,1,4:5), set(gca,'fontsize',12), hold on
% plot(Xord(iobs)*kpd,[eqar_gd(iobs).dT],'.g','MarkerSize',10)
% plot(Xord(iobs)*kpd,[eqar_gd(iobs).dT_comb],'.b','MarkerSize',18)
% plot(([eqar_gd(ilan).slon]+124)*kpd + 440,[eqar_gd(ilan).dT],'.g','MarkerSize',10)
% plot(([eqar_gd(ilan).slon]+124)*kpd + 440,[eqar_gd(ilan).dT_comb],'.r','MarkerSize',18)
% % plot(Xord(ilan)*kpd,[eqar_gd(ilan).dT],'.r','MarkerSize',18)
% % figure things
% grid on
% text(-115,0,'$\Delta t_S$','Fontsize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex','Rotation',90)
% text(320,-4.5,'Distance from ridge (km)','Fontsize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex')
% xlim([-80 720]), ylim([-3 3])

%% ===================================================================
%% ---------------------------- WAVEFORMS ----------------------------
%% ===================================================================

for is = 1:length(eqar_gd)
    dat = eqar_gd(is).(['dat',comp])';

    W = 2*filtfs./eqar_gd(is).samprate;
    [fb, fa] = butter(2,W);
    dd(:,is) = filtfilt(fb,fa,dat);
end

%% ------- plot normalised to one value
figure(2), clf, set(gcf,'position',[150 000 600 800]), hold on
for is = 1:length(eqar_gd)
	tt = (eqar_gd(is).tt-eqar_gd(is).abs_arrT)*spd;

    hp = plot(tt,5*dd(:,is)/max(max(abs(dd))) + Xord(is),'LineWidth',2);
    set(hp,'color',colour_get(abs(eqar_gd(is).Xrft),Xrftlims_plot(2),Xrftlims_plot(1),flipud(jet)))
    text(datwind(1)+2,Xord(is)+0.4,eqar_gd(is).sta,...
        'FontSize',8,'interpreter','latex','HorizontalAlignment','left')
    text(datwind(2) + 1,Xord(is),sprintf('%.0f km',eqar_gd(is).Xrft*kpd),...
        'FontSize',10,'interpreter','latex')
end
plot([1;1]*window,[0;length(eqar_gd)+1]*[1 1],'--k','LineWidth',1)
set(gca,'FontSize',12,'YTick',[],'position',[0.16 0.11 0.75 0.815])
axis([datwind,-0.5,length(eqar_gd)+1.5])
xlabel('Time from phase onset (s)','FontSize',14,'interpreter','latex')
title(sprintf('%s $~$ %s-wave, %s-comp, $f_{hi}~$: %.2f',datestr(evtime,'yyyy-mm-dd HH:MM:SS'),...
    phase,comp,fmax),'Fontsize',16,'Interpreter','latex')

%% -------  plot self normalized
figure(22), clf, set(gcf,'position',[750 000 600 800]), hold on
for is = 1:length(eqar_gd)
	tt = (eqar_gd(is).tt-eqar_gd(is).abs_arrT)*spd;

    hp = plot(tt,2*dd(:,is)/max(abs(dd(:,is))) + Xord(is),'LineWidth',2);
    set(hp,'color',colour_get(abs(eqar_gd(is).Xrft),Xrftlims_plot(2),Xrftlims_plot(1),flipud(jet)))
    text(datwind(1)+2,Xord(is)+0.4,eqar_gd(is).sta,...
        'FontSize',8,'interpreter','latex','HorizontalAlignment','left')
    text(datwind(2) + 1,Xord(is),sprintf('%.0f km',eqar_gd(is).Xrft*kpd),...
        'FontSize',10,'interpreter','latex')
end
plot([1;1]*window,[0;length(eqar_gd)+1]*[1 1],'--k','LineWidth',1)
set(gca,'FontSize',12,'YTick',[],'position',[0.16 0.11 0.75 0.815])
axis([datwind,-0.5,length(eqar_gd)+1.5])
xlabel('Time from phase onset (s)','FontSize',14,'interpreter','latex')
title(sprintf('%s $~$ %s-wave, %s-comp, $f_{hi}~$: %.2f',datestr(evtime,'yyyy-mm-dd HH:MM:SS'),...
    phase,comp,fmax),'Fontsize',16,'Interpreter','latex')

% % ridge axis label
% text(datwind(1)-2,sum([eqar_gd.Xord] < 0)+0.5,'\textbf{Axis $\succ$}',...
%     'HorizontalAlignment','right','VerticalAlignment','bottom',...
%     'FontSize',15,'FontWeight','bold','Interpreter','latex')
% % deformation front label
% text(datwind(1)-2,sum([eqar_gd.Xord] < 305/kpd)+0.5,'\textbf{DF $\succ$}',...
%     'HorizontalAlignment','right','VerticalAlignment','bottom',...
%     'FontSize',15,'FontWeight','bold','Interpreter','latex')
% % Coastline label
% text(datwind(1)-2,sum([eqar_gd.Xord] < 440/kpd)+0.5,'\textbf{Coast $\succ$}',...
%     'HorizontalAlignment','right','VerticalAlignment','bottom',...
%     'FontSize',15,'FontWeight','bold','Interpreter','latex')

%% ===================================================================
%% ------------------------ SPECTRA (W/ FITS) ------------------------
%% ===================================================================
figure(3), clf, set(gcf,'position',[200 10 800 600]), hold on
for is = 1:length(eqar_gd)
    
    frq = eqar_gd(is).frq;
    lnR = log(eqar_gd(is).specss./specss_ref);

    ind = frq <= eqar_gd(is).fcrosshi & frq >= eqar_gd(is).fcrosslo;

    fo = fit(frq(ind),lnR(ind),'poly1');
  
    hp = plot(frq(ind),lnR(ind) - fo.p2,'o-','LineWidth',1.5);
    hpf = plot(frq(ind),fo.p1*frq(ind),'--');

    set([hp,hpf],'color',colour_get(abs(eqar_gd(is).Xrft),Xrftlims_plot(2),Xrftlims_plot(1),flipud(jet)))
%     xlim([0 par.hifrq])
    
end
xlim([min([eqar_gd.fcrosslo])-0.1 max([eqar_gd.fcrosshi])+0.1])
% ylim([-2.5 3])
set(gca,'FontSize',8)
xlabel('Frequency (Hz)','interpreter','latex','FontSize',16)
ylabel('Spectral ratio: $~\ln~(A_i/A_0)$','interpreter','latex','FontSize',16)
title('OBS spectral ratios relative to mean spectrum','interpreter','latex','FontSize',18)
legend off

cbar_custom(gca, 'location',[0.06 .12 1 1.18],'tickside','bottom',...
    'lims',kpd*Xrftlims_plot,'tickvals',round_level([kpd*Xrftlims_plot(1):50:kpd*Xrftlims_plot(2)],100),...
    'FontSize',10,'FontWeight','bold','cmap',flipud(jet),...
    'title','Distance E of rift (km)','interpreter','latex');




end

