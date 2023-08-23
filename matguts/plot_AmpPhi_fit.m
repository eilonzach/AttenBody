function [misfit_amp,misfit_phi,E_a,E_p] = plot_AmpPhi_fit( Amp,Phi,fmids,wts,dtstar,dT,A0, alp,ghandl,w0 )
% [misfit_amp,misfit_phi] = plot_AmpPhi_fit( Amp,Phi,fmids,[wts=1],dtstar,dT,A0, [alpha=0], [ghandl=33], [w0=2pi] )
% 
% Quick function to plot amplitude and phase spectra and the fits to them,
% along with the misfit.

if nargin < 8 || isempty(alp)
    alp = 0;
end
if nargin < 9 || isempty(ghandl)
    ghandl = 33;
end
if nargin < 10 || isempty(w0)
    w0 = 2*pi;
end
if isempty(wts)
    wts = ones(size(Amp));
end

Amp = Amp(:);
Phi = Phi(:);
fmids = fmids(:);
ws = 2*pi*fmids;
f0 = w0/2/pi;

ffa = fmids;
ffp = fmids;
xlaba = 'freq (Hz)';
xlabq = 'freq (Hz)';

ffall = linspace(0.033,1,100);

if alp == 0
    lnApred   = log(A0) - pi*ffa*dtstar;
    phipred = -(1/pi)*log(ffa)*dtstar + dtstar + dT;
%     ffa = fmids;
%     ffp = fmids;
%     xlaba = 'freq (Hz)';
%     xlabq = 'freq (Hz)';
    lnApredall   = log(A0) - pi*ffall*dtstar;
    phipredall = -(1/pi)*log(ffall)*dtstar + dtstar + dT;
else
    lnApred = log(A0) - 0.5*ws.^(1-alp)*w0^(alp)*dtstar;
    phipred = 0.5*(ws/w0).^(-alp) * cot(alp*pi/2)*dtstar + dT;
%     ffa = fmids.^(1-alp);
%     ffp = fmids.^(-alp);
%     xlaba = 'freq$^{(1-\alpha)}$';
%     xlabq = 'freq$^{(-\alpha)}$';
    lnApredall = log(A0) - pi*(ffall).^(1-alp)*f0^(alp)*dtstar;
    phipredall = 0.5*(ffall/f0).^(-alp) * cot(alp*pi/2)*dtstar + dT;
end

E_a = Amp-exp(lnApred);
E_p = Phi-phipred;

misfit_amp = E_a'*diag(wts)*E_a;
misfit_phi = E_p'*diag(wts)*E_p;


%% START THE FANCY PLOTS
if length(ghandl) == 1
    fig = figure(ghandl); clf(fig), set(fig,'pos',[600 600 600,800])
    ax1 = subplot(2,1,1);
    ax2 = subplot(2,1,2);
elseif length(ghandl) == 2
    ax1 = ghandl(1);
    ax2 = ghandl(2);
end

wts(wts==0) = nan;

% ================ PLOT THE AMPLITUDE SPECTRA ================
axes(ax1), hold on

h = scatter(ax1,fmids,log(Amp),150*wts); set(h,'MarkeredgeColor','k','MarkerFaceColor','r','linewidth',1.5)
plot(ax1,ffall,lnApredall,':c','Linewidth',1.5)
plot(ax1,ffa,lnApred,'g','Linewidth',2.5)

xlabel(ax1,xlaba,'FontSize',22,'interpreter','latex')
ylabel(ax1,'$\ln\,(R_{12})$','FontSize',22,'interpreter','latex')
title(ax1,sprintf('$\\Delta t^*$ = %.2f \\,\\,\\, $\\alpha$ = %.2f',dtstar,alp),'FontSize',22,'interpreter','latex')
set(ax1,'fontsize',15,'Xscale','linear','ylim',[-3 1],'linewidth',2,'box','on','xlim',[0 1])

ax = axis; xl = ax(1) + 0.05*(ax(2)-ax(1)); yl = ax(4) - 0.1*(ax(4)-ax(3));
text(ax1,xl,yl,['Weighted misfit = ',num2str(misfit_amp)],'interpreter','latex','fontsize',20,'horizontalalignment','left','verticalalignment','bottom')


% ================ PLOT THE PHASE SPECTRA ================
axes(ax2), hold on
scatter(ax2,fmids,Phi,150*wts,'ok','MarkerFaceColor','r','linewidth',1.5)
plot(ax2,ffall,phipredall,':c','Linewidth',1.5)
plot(ax2,ffp,phipred,'g','Linewidth',2.5)

% ha1 = arrow([eqar(ind1).fcross,1],[eqar(ind1).fcross,0],5,90,2,0.1,'FaceColor','m'); % fcross for sta 1
% ha2 = arrow([eqar(ind2).fcross,1],[eqar(ind2).fcross,0],5,90,2,0.1,'FaceColor','g'); % fcross for sta 2
% plot(fmax*[1 1],[-1 1],'--b')  

xlabel(ax2,xlabq,'FontSize',22,'interpreter','latex')
ylabel(ax2,'$\Delta \phi_{12}$ (s)','FontSize',22,'interpreter','latex')
set(ax2,'fontsize',15,'Xscale','linear','linewidth',2,'box','on','xlim',[0 1])

ax = axis; xl = ax(1) + 0.05*(ax(2)-ax(1)); yl = ax(4) - 0.1*(ax(4)-ax(3));
text(xl,yl,['Weighted misfit = ',num2str(misfit_phi)],'interpreter','latex','fontsize',20,'horizontalalignment','left','verticalalignment','bottom')

return
% ================ SAVE FIGURE ================
ofile = sprintf('synth_comb_Q1-%.0f_Q2-%.0f',Q0_1,Q0_2);
save2pdf(33,ofile,'figs');


end

