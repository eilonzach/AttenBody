function plot_trans_func(method_trans_func,trans_func)
%PLOT_TRANS_FUNC Summary of this function goes here
%   Detailed explanation goes here

fs = fieldnames(trans_func);
for ii = 1:length(fs), eval(sprintf('%s = trans_func.%s;',fs{ii},fs{ii})); end

pper = 1/16;
sper = 1/8;

figure(111), clf, set(gcf,'pos',[10 10 1200 790])

switch method_trans_func
    case 1
        subplot(611), hold on
        loglog(f,smooth(abs(lpz_12),40),'b-');
        loglog(f,smooth(abs(lpz_12_use),40),'g-');
        title([sta,' Transfer Function, lpz-12 ']);
        plot([pper pper],[min(abs(lpz_12))*10 max(abs(lpz_12))/10],'-b','LineWidth',2)
        plot([sper sper],[min(abs(lpz_12))*10 max(abs(lpz_12))/10],'r','LineWidth',2)
        plot([freqcomp freqcomp],[min(abs(lpz_12))*10 max(abs(lpz_12))/10],'-k')
        set(gca,'ylim',[min(abs(lpz_12))/2 max(abs(lpz_12))*2],'xscale','log','yscale','log')

        subplot(612), hold on
        loglog(f,smooth(abs(l1z),40),'b-');
        loglog(f,smooth(abs(l1z_use),40),'g-');
        title([sta,' Transfer Function, l1z ']);
        plot([pper pper],[min(abs(l1z))*10 max(abs(l1z))/10],'-b','LineWidth',2)
        plot([sper sper],[min(abs(l1z))*10 max(abs(l1z))/10],'r','LineWidth',2)
        plot([freqcomp freqcomp],[min(abs(l1z))*10 max(abs(l1z))/10],'-k')
        set(gca,'ylim',[min(abs(l1z))/2 max(abs(l1z))*2],'xscale','log','yscale','log')

        subplot(613), hold on
        loglog(f,smooth(abs(l2z_1),40),'b-');
        loglog(f,smooth(abs(l2z_1_use),40),'g-');
        title([sta,' Transfer Function, l2z-1 ']);
        plot([pper pper],[min(abs(l2z_1))*10 max(abs(l2z_1))/10],'-b','LineWidth',2)
        plot([sper sper],[min(abs(l2z_1))*10 max(abs(l2z_1))/10],'r','LineWidth',2)
        plot([freqcomp freqcomp],[min(abs(l2z_1))*10 max(abs(l2z_1))/10],'-k')
        set(gca,'ylim',[min(abs(l2z_1))/2 max(abs(l2z_1))*2],'xscale','log','yscale','log')

        subplot(614), hold on
        loglog(f,smooth(abs(l12),40),'b-');
        loglog(f,smooth(abs(l12_use),40),'g-');
        title([sta,' Transfer Function, l12 ']);
        plot([pper pper],[min(abs(l12))*10 max(abs(l12))/10],'-b','LineWidth',2)
        plot([sper sper],[min(abs(l12))*10 max(abs(l12))/10],'r','LineWidth',2)
        plot([freqcomp freqcomp],[min(abs(l12))*10 max(abs(l12))/10],'-k')
        set(gca,'ylim',[min(abs(l12))/2 max(abs(l12))*2],'xscale','log','yscale','log')

        subplot(615), hold on
        loglog(f,smooth(abs(l1p),40),'b-');
        loglog(f,smooth(abs(l1p_use),40),'g-');
        title([sta,' Transfer Function, l1p ']);
        plot([pper pper],[min(abs(l1p))*10 max(abs(l1p))/10],'-b','LineWidth',2)
        plot([sper sper],[min(abs(l1p))*10 max(abs(l1p))/10],'r','LineWidth',2)
        plot([freqcomp freqcomp],[min(abs(l1p))*10 max(abs(l1p))/10],'-k')
        set(gca,'ylim',[min(abs(l1p))/2 max(abs(l1p))*2],'xscale','log','yscale','log')

        subplot(616), hold on
        loglog(f,smooth(abs(l2p_1),40),'b-');
        loglog(f,smooth(abs(l2p_1_use),40),'g-');
        title([sta,' Transfer Function, l2p-1 ']);
        plot([pper pper],[min(abs(l2p_1))*10 max(abs(l2p_1))/10],'-b','LineWidth',2)
        plot([sper sper],[min(abs(l2p_1))*10 max(abs(l2p_1))/10],'r','LineWidth',2)
        plot([freqcomp freqcomp],[min(abs(l2p_1))*10 max(abs(l2p_1))/10],'-k')
        set(gca,'ylim',[min(abs(l2p_1))/2 max(abs(l2p_1))*2],'xscale','log','yscale','log')
        
    case 2 
        subplot(311), hold on
        loglog(f,smooth(abs(lpz),40),'b-');
        loglog(f,smooth(abs(lpz_use),40),'g-');
        title([sta,' Transfer Function, lpz ']);
        plot([pper pper],[min(abs(lpz))*10 max(abs(lpz))/10],'-b','LineWidth',2)
        plot([sper sper],[min(abs(lpz))*10 max(abs(lpz))/10],'r','LineWidth',2)
        plot([freqcomp freqcomp],[min(abs(lpz))*10 max(abs(lpz))/10],'-k')
        set(gca,'ylim',[min(abs(lpz))/2 max(abs(lpz))*2],'xscale','log','yscale','log')

        subplot(312), hold on
        loglog(f,smooth(abs(lp1),40),'b-');
        title([sta,' Transfer Function, lp1 ']);    
        plot([pper pper],[min(abs(lp1))*10 max(abs(lp1))/10],'-b','LineWidth',2)
        plot([sper sper],[min(abs(lp1))*10 max(abs(lp1))/10],'r','LineWidth',2)
        plot([freqcomp freqcomp],[min(abs(lp1))*10 max(abs(lp1))/10],'-k')
        set(gca,'ylim',[min(abs(lp1))/2 max(abs(lp1))*2],'xscale','log','yscale','log')

        subplot(313), hold on
        loglog(f,smooth(abs(lp2),40),'b-');
        title([sta,' Transfer Function, lp2 ']);
        plot([pper pper],[min(abs(lp2))*10 max(abs(lp2))/10],'-b','LineWidth',2)
        plot([sper sper],[min(abs(lp2))*10 max(abs(lp2))/10],'r','LineWidth',2)
        plot([freqcomp freqcomp],[min(abs(lp2))*10 max(abs(lp2))/10],'-k')
        set(gca,'ylim',[min(abs(lp2))/2 max(abs(lp2))*2],'xscale','log','yscale','log')
    
    case 3
        subplot(611)
        loglog(f,smooth(abs(l1z_p2),40),'b-');
        title([sta,' Transfer Function, l1z-p2 ']);
        hold on
        plot([pper pper],[min(abs(l1z_p2))*10 max(abs(l1z_p2))/10],'-b','LineWidth',2)
        hold on
        plot([sper sper],[min(abs(l1z_p2))*10 max(abs(l1z_p2))/10],'r','LineWidth',2)
        hold on
        plot([freqcomp freqcomp],[min(abs(l1z_p2))*10 max(abs(l1z_p2))/10],'-k')
        hold on
        subplot(612)
        loglog(f,smooth(abs(lpz),40),'b-');
        title([sta,' Transfer Function, lpz ']);
        hold on
        plot([pper pper],[min(abs(lpz))*10 max(abs(lpz))/10],'-b','LineWidth',2)
        hold on
        plot([sper sper],[min(abs(lpz))*10 max(abs(lpz))/10],'r','LineWidth',2)
        hold on
        plot([freqcomp freqcomp],[min(abs(lpz))*10 max(abs(lpz))/10],'-k')
        hold on
        subplot(613)
        loglog(f,smooth(abs(l2z_p),40),'b-');
        title([sta,' Transfer Function, l2z-p ']);
        hold on
        plot([pper pper],[min(abs(l2z_p))*10 max(abs(l2z_p))/10],'-b','LineWidth',2)
        hold on
        plot([sper sper],[min(abs(l2z_p))*10 max(abs(l2z_p))/10],'r','LineWidth',2)
        hold on
        plot([freqcomp freqcomp],[min(abs(l2z_p))*10 max(abs(l2z_p))/10],'-k')
        hold on
        subplot(614)
        loglog(f,smooth(abs(lp2),40),'b-');
        title([sta,' Transfer Function, lp2 ']);
        hold on
        plot([pper pper],[min(abs(lp2))*10 max(abs(lp2))/10],'-b','LineWidth',2)
        hold on
        plot([sper sper],[min(abs(lp2))*10 max(abs(lp2))/10],'r','LineWidth',2)
        hold on
        plot([freqcomp freqcomp],[min(abs(lp2))*10 max(abs(lp2))/10],'-k')
        hold on
        subplot(615)
        loglog(f,smooth(abs(lp1),40),'b-');
        title([sta,' Transfer Function, lp1 ']);
        hold on
        plot([pper pper],[min(abs(lp1))*10 max(abs(lp1))/10],'-b','LineWidth',2)
        hold on
        plot([sper sper],[min(abs(lp1))*10 max(abs(lp1))/10],'r','LineWidth',2)
        hold on
        plot([freqcomp freqcomp],[min(abs(lp1))*10 max(abs(lp1))/10],'-k')
        hold on
        subplot(616)
        loglog(f,smooth(abs(l21_p),40),'b-');
        title([sta,' Transfer Function, l21-p ']);
        hold on
        plot([pper pper],[min(abs(l21_p))*10 max(abs(l21_p))/10],'-b','LineWidth',2)
        hold on
        plot([sper sper],[min(abs(l21_p))*10 max(abs(l21_p))/10],'r','LineWidth',2)
        hold on
        plot([freqcomp freqcomp],[min(abs(l21_p))*10 max(abs(l21_p))/10],'-k')
        hold on
    case 4
        fprintf('Sorry, this method not yet ready for primetime!\n'), return
end

% %% SAVE FIGURES
% figure(101)
% set(gcf,'PaperPositionMode','manual');
% set(gcf,'PaperUnits','inches');
% set(gcf,'PaperOrientation','portrait');
% set(gcf,'PaperPosition',[.05 .05 8 10.5]);
% if method ==1
%     filename=sprintf('%s/TransferFunctions_met1corrected_%s',figoutpath,station);
% elseif method ==2
%     filename=sprintf('%s/TransferFunctions_met2corrected_%s',figoutpath,station);
%     elseif method ==3
%     filename=sprintf('%s/TransferFunctions_met3corrected_%s',figoutpath,station);
% end
% print(gcf,'-dpng',filename)

end

