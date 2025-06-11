function figh = plot_focustest_pretty(stas,Amat,phimat,fmids,wtmat,test_a2pwt,a2ptest_dtstar,a2ptest_A0,a2ptest_Eamp_stawise,a2ptest_Ephi_stawise,focusscore,focus_thresh)


%% some prelims
Nstas = length(stas);
Nf    = length(fmids);

% set up figure
figh = figure(32); clf,set(figh,'position',[1406 148 985 1179],'color','w');
% axes
ax_lnA = axes('position',[0.0609 0.7006 0.4162 0.2480],'layer','top','box','on'); hold on
ax_phi = axes('position',[0.0609 0.3706 0.4162 0.2480],'layer','top','box','on'); hold on
ax_dts = axes('position',[0.5609 0.7006 0.4162 0.2480],'layer','top','box','on','xscale','log'); hold on
ax_dA0 = axes('position',[0.5609 0.3706 0.4162 0.2480],'layer','top','box','on','xscale','log'); hold on
ax_Err = axes('position',[0.5609 0.0506 0.4162 0.2480],'layer','top','box','on','xscale','log'); hold on
ax_foc = axes('position',[0.0609 0.0506 0.4162 0.2480],'layer','top','box','on'); hold on

% colours
fs_use = focusscore+1; fs_use(isinf(fs_use)) = 999;

stacol = colour_get(log(fs_use),ceil(max(log(fs_use))),floor(min(log(fs_use))),[brewermap(20,'PuBuGn')]).^1.5;
% stacol = colour_get([1:Nstas]',Nstas,1,[brewermap(10,'OrRd');brewermap(10,'RdPu');brewermap(20,'PuBuGn')]).^1.5;
% stacol = colour_get(a2ptest_dtstar(:,test_a2pwt==2),1.5,-1.5,[brewermap(10,'OrRd');brewermap(10,'RdPu');brewermap(20,'PuBuGn')]).^1.5;
% stacol = brewermap(Nstas,'Dark2');

%% compute spectra relative to mean spectrum - have to invert the pairwise. 
% Good news is each frequency is independent. loop through them
% first make G matrix for handshake
Npair = size(Amat,1);
G = zeros(Npair+1,Nstas);
count = 0;
for is1 = 1:Nstas   
for is2 = is1+1:Nstas
    count = count+1;
    G(count,[is1,is2]) = [-1 1]; % hope it's this way round.... Will see at the end.
end
end
G(end,:) = 1; % constraint that average differential specturum is zero
% now build single-station (i.e. mean-ed) differential spectra
lnA_av_f = nan(Nstas,Nf);
phi_av_f = nan(Nstas,Nf);
wt_av_f = nan(Nstas,Nf);
for iff = 1:Nf
    inoignore = wtmat(:,iff)~=0; % don't ignore if non-zero weight (or constraint!)
    G_ = G([inoignore;true],:);
    Ginv = (G_'*G_ + 1e-9*eye(Nstas))\G_';
    lnA_av_f(:,iff) = Ginv*[log(Amat(inoignore,iff));0]; % m = Finv*[d;h]
    phi_av_f(:,iff) = Ginv*[phimat(inoignore,iff);0]; % m = Finv*[d;h]
    wt_av_f(:,iff) = Ginv*[wtmat(inoignore,iff);sum(wtmat(inoignore,iff))]; % m = Finv*[d;h]
end


% while testing
% count = 1;
% for is1 = 1:Nstas
%     for is2 = is1+1:Nstas
%         plot(ax_lnA,fmids,Amat(count,:),'color',stacol(is1,:),'linewidth',1)
%         plot(ax_phi,fmids,phimat(count,:),'color',stacol(is1,:),'linewidth',1)
%         plot(ax_lnA,fmids,Amat(count,:)-.01,'color',stacol(is2,:))
%         plot(ax_phi,fmids,phimat(count,:)-.01,'color',stacol(is2,:))
%         count = count+1;
%     end 
%     plot(ax_dts,fmids,exp(lnA_av_f(is1,:)),'color',stacol(is1,:))
%     plot(ax_dA0,fmids,phi_av_f(is1,:),'color',stacol(is1,:))
% end



%% plot spectra
for is1 = 1:Nstas
    plot(ax_lnA,fmids,exp(lnA_av_f(is1,:)),'color',stacol(is1,:),'LineWidth',2)
    plot(ax_phi,fmids,phi_av_f(is1,:),'color',stacol(is1,:),'LineWidth',2)
    try scatter(ax_lnA,fmids,exp(lnA_av_f(is1,:)),60*wt_av_f(is1,:)+1e-9,stacol(is1,:),'filled','markeredgecolor','k'); 
    catch, continue; end
    try scatter(ax_phi,fmids,phi_av_f(is1,:),60*wt_av_f(is1,:)+1e-9,stacol(is1,:),'filled','markeredgecolor','k'); catch, continue; end
end

%% plot outcome of A2phiwt tests
for is1 = 1:Nstas
    plot(ax_dts,test_a2pwt,a2ptest_dtstar(is1,:),'color',stacol(is1,:),'LineWidth',2)
    plot(ax_dA0,test_a2pwt,log(a2ptest_A0(is1,:)),'color',stacol(is1,:),'LineWidth',2)
    plot(ax_Err,test_a2pwt,a2ptest_Eamp_stawise(is1,:),'-','color',stacol(is1,:),'LineWidth',2)
    plot(ax_Err,test_a2pwt,a2ptest_Ephi_stawise(is1,:),'--','color',stacol(is1,:),'LineWidth',2)
end

%% plot focus score
[~,isort] = sort(focusscore);
for isloop = 1:Nstas
    is1 = isort(isloop);
    plot(ax_foc,isloop,focusscore(is1),'p','markersize',20,'MarkerfaceColor',stacol(is1,:),'MarkeredgeColor',stacol(is1,:).^3,'LineWidth',1.5)
end
yline(focus_thresh,'--r')

%% axes prettying/labelling
set([ax_lnA,ax_phi,ax_dts,ax_dA0,ax_Err,ax_foc],'linewidth',1.8,'fontsize',15)
set([ax_lnA,ax_phi],'xlim',[0.045 0.5])
% special
set(ax_phi,'xscale','log','xtick',[0.05,1/10,1/5,1/2,1])
set(ax_foc,'xtick',1:Nstas,'XTickLabel',stas(isort));

% labels
xlabel(ax_lnA,'frequency (Hz)','FontSize',18,'FontWeight','bold')
xlabel(ax_phi,'log(frequency (Hz))','FontSize',18,'FontWeight','bold')
xlabel(ax_foc,'station','FontSize',18,'FontWeight','bold')
xlabel([ax_dts,ax_dA0,ax_Err],'Amp2Phi weight','FontSize',18,'FontWeight','bold')
ylabel(ax_lnA,'Relative amplitude','FontSize',18,'FontWeight','bold')
ylabel(ax_phi,'Differential phase (s)','FontSize',18,'FontWeight','bold')
ylabel(ax_foc,'Focus score','FontSize',18,'FontWeight','bold')
ylabel(ax_dts,'\delta t^*','FontSize',18,'FontWeight','bold')
ylabel(ax_dA0,'\delta A_0','FontSize',18,'FontWeight','bold')
ylabel(ax_Err,'Error in Amp (—) or phase (--)','FontSize',18,'FontWeight','bold')




end