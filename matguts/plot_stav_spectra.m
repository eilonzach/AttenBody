% code run after evdata_4a_STAV_SPECTRA that plots average spectral properties for a single sta

% inherited from outside:
%       all stacks
%       mean and median values
%       ifsavefigs option
%       odir

odir = [specdir,sta,'/'];

if ifsavefigs
    if exist(odir,'dir')~=7, mkdir(odir); end
end

grey = [0.5 0.5 0.5];

pper = 1/16;
sper = 1/8;

T = (win_pt_end(1)-win_pt_start(1)+1)*dt;

%% plots of averages

%% Power plots
figure(3), clf
subplot(411); hold on
title(sprintf('Z-component, Station: %s',sta));
for ik = 1:kk
    loglog(f,smooth(czz_all(ik,:),40),'-','LineWidth',.5,'Color',grey);
end
loglog(f,smooth(czz_mean,40),'-','LineWidth',.7,'Color','g');
loglog(f,smooth(czz_med,40),'-','LineWidth',.7,'Color','c');
xlabel('Frequency (Hz)')
ylabel('Power Spectrum')
% ylim([10^-20 max(maxpowz)])
xlim([2./T .5./dt])
plot([pper pper],max(maxpowz)*[1e-10 1],'-b','LineWidth',2)
plot([sper sper],max(maxpowz)*[1e-10 1],'r','LineWidth',2)
plot([freqcomp freqcomp],max(maxpowz)*[1e-10 1],'-k','LineWidth',2)
set(gca,'XScale','log','Yscale','log')

subplot(412); hold on
title(sprintf('H1-component, Station: %s',sta));
for ik = 1:kk
    loglog(f,smooth(c11_all(ik,:),40),'-','LineWidth',.5,'Color',grey);
end
loglog(f,smooth(c11_mean,40),'-','LineWidth',.7,'Color','g');
loglog(f,smooth(c11_med,40),'-','LineWidth',.7,'Color','c');
xlabel('Frequency (Hz)')
ylabel('Power Spectrum')
% ylim([10^-20 max(maxpowz)])
xlim([2./T .5./dt])
plot([pper pper],max(maxpow1)*[1e-10 1],'-b','LineWidth',2)
plot([sper sper],max(maxpow1)*[1e-10 1],'r','LineWidth',2)
plot([freqcomp freqcomp],max(maxpow1)*[1e-10 1],'-k','LineWidth',2)
set(gca,'XScale','log','Yscale','log')

subplot(413); hold on
title(sprintf('H2-component, Station: %s',sta));
for ik = 1:kk
    loglog(f,smooth(c22_all(ik,:),40),'-','LineWidth',.5,'Color',grey);
end
loglog(f,smooth(c22_mean,40),'-','LineWidth',.7,'Color','g');
loglog(f,smooth(c22_med,40),'-','LineWidth',.7,'Color','c');
xlabel('Frequency (Hz)')
ylabel('Power Spectrum')
% ylim([10^-20 max(maxpowz)])
xlim([2./T .5./dt])
plot([pper pper],max(maxpow2)*[1e-10 1],'-b','LineWidth',2)
plot([sper sper],max(maxpow2)*[1e-10 1],'r','LineWidth',2)
plot([freqcomp freqcomp],max(maxpow2)*[1e-10 1],'-k','LineWidth',2)
set(gca,'XScale','log','Yscale','log')

subplot(414); hold on
title(sprintf('P-component, Station: %s',sta));
for ik = 1:kk
    loglog(f,smooth(cpp_all(ik,:),40),'-','LineWidth',.5,'Color',grey);
end
loglog(f,smooth(cpp_mean,40),'-','LineWidth',.7,'Color','g');
loglog(f,smooth(cpp_med,40),'-','LineWidth',.7,'Color','c');
xlabel('Frequency (Hz)')
ylabel('Power Spectrum')
% ylim([10^-20 max(maxpowz)])
xlim([2./T .5./dt])
plot([pper pper],max(maxpowp)*[1e-10 1],'-b','LineWidth',2)
plot([sper sper],max(maxpowp)*[1e-10 1],'r','LineWidth',2)
plot([freqcomp freqcomp],max(maxpowp)*[1e-10 1],'-k','LineWidth',2)
set(gca,'XScale','log','Yscale','log')

set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);
if ifsavefigs
    filename=[odir,sta,'_STAV_power'];
    print(gcf,'-dpng',filename)
end 

%% coherence
figure(4), clf
subplot(611); hold on
for ik = 1:kk
    semilogx(f,smooth(coh1z_all(ik,:),40),'-','LineWidth',.5,'Color',grey);
end
semilogx(f,smooth(coh1z_mean,40),'-','LineWidth',.7,'Color','g');
title(sprintf('Z-1, Station: %s',sta))
xlabel('Frequency (Hz)')
ylabel('Coherence')
ylim([0 1])
xlim([2./T .5./dt])
semilogx([pper pper],[0 1],'-b','LineWidth',2)
semilogx([sper sper],[0 1],'-r','LineWidth',2)
semilogx([freqcomp freqcomp],[0,1],'-k','LineWidth',2)
set(gca,'XScale','log')

subplot(612); hold on
for ik = 1:kk
    semilogx(f,smooth(coh2z_all(ik,:),40),'-','LineWidth',.5,'Color',grey);
end
semilogx(f,smooth(coh2z_mean,40),'-','LineWidth',.7,'Color','g');
title(sprintf('Z-2, Station: %s',sta))
xlabel('Frequency (Hz)')
ylabel('Coherence')
ylim([0 1])
xlim([2./T .5./dt])
semilogx([pper pper],[0 1],'-b','LineWidth',2)
semilogx([sper sper],[0 1],'-r','LineWidth',2)
semilogx([freqcomp freqcomp],[0,1],'-k','LineWidth',2)
set(gca,'XScale','log')

subplot(613); hold on
for ik = 1:kk
    semilogx(f,smooth(cohpz_all(ik,:),40),'-','LineWidth',.5,'Color',grey);
end
semilogx(f,smooth(cohpz_mean,40),'-','LineWidth',.7,'Color','g');
title(sprintf('P-Z, Station: %s',sta))
xlabel('Frequency (Hz)')
ylabel('Coherence')
ylim([0 1])
xlim([2./T .5./dt])
semilogx([pper pper],[0 1],'-b','LineWidth',2)
semilogx([sper sper],[0 1],'-r','LineWidth',2)
semilogx([freqcomp freqcomp],[0,1],'-k','LineWidth',2)
set(gca,'XScale','log')

subplot(614); hold on
for ik = 1:kk
    semilogx(f,smooth(coh12_all(ik,:),40),'-','LineWidth',.5,'Color',grey);
end
semilogx(f,smooth(coh12_mean,40),'-','LineWidth',.7,'Color','g');
title(sprintf('1-2, Station: %s',sta))
xlabel('Frequency (Hz)')
ylabel('Coherence')
ylim([0 1])
xlim([2./T .5./dt])
semilogx([pper pper],[0 1],'-b','LineWidth',2)
semilogx([sper sper],[0 1],'-r','LineWidth',2)
semilogx([freqcomp freqcomp],[0,1],'-k','LineWidth',2)
set(gca,'XScale','log')

subplot(615); hold on
for ik = 1:kk
    semilogx(f,smooth(coh1p_all(ik,:),40),'-','LineWidth',.5,'Color',grey);
end
semilogx(f,smooth(coh1p_mean,40),'-','LineWidth',.7,'Color','g');
title(sprintf('1-P, Station: %s',sta))
xlabel('Frequency (Hz)')
ylabel('Coherence')
ylim([0 1])
xlim([2./T .5./dt])
semilogx([pper pper],[0 1],'-b','LineWidth',2)
semilogx([sper sper],[0 1],'-r','LineWidth',2)
semilogx([freqcomp freqcomp],[0,1],'-k','LineWidth',2)
set(gca,'XScale','log')

subplot(616); hold on
for ik = 1:kk
    semilogx(f,smooth(coh2p_all(ik,:),40),'-','LineWidth',.5,'Color',grey);
end
semilogx(f,smooth(coh2p_mean,40),'-','LineWidth',.7,'Color','g');
title(sprintf('2-P, Station: %s',sta))
xlabel('Frequency (Hz)')
ylabel('Coherence')
ylim([0 1])
xlim([2./T .5./dt])
semilogx([pper pper],[0 1],'-b','LineWidth',2)
semilogx([sper sper],[0 1],'-r','LineWidth',2)
semilogx([freqcomp freqcomp],[0,1],'-k','LineWidth',2)
set(gca,'XScale','log')

set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);
if ifsavefigs
    filename=[odir,sta,'_STAV_coherence'];
    print(gcf,'-dpng',filename)
end

%     figure(6)
%     subplot(611);
%     title(sprintf('Z-1, Station: %s',sta))
%     xlabel('Frequency (Hz)')
%     ylabel('Phase')
%     ylim([-90 90])
%     plot([pper pper],[-90 90],'-b','LineWidth',2)
%     hold on
%     plot([sper sper],[-90 90],'-r','LineWidth',2)
%     hold on
%     plot([freqcomp freqcomp],[-90 90],'-k','LineWidth',2)
%     hold on
%     subplot(612);
%     title(sprintf('Z-2, Station: %s',sta))
%     xlabel('Frequency (Hz)')
%     ylabel('Phase')
%     ylim([-90 90])
%     plot([pper pper],[-90 90],'-b','LineWidth',2)
%     hold on
%     plot([sper sper],[-90 90],'-r','LineWidth',2)
%     hold on
%     plot([freqcomp freqcomp],[-90 90],'-k','LineWidth',2)
%     hold on
%     subplot(613);
%     title(sprintf('P-Z, Station: %s',sta))
%     xlabel('Frequency (Hz)')
%     ylabel('Phase')
%     ylim([-90 90])
%     plot([pper pper],[-90 90],'-b','LineWidth',2)
%     hold on
%     plot([sper sper],[-90 90],'-r','LineWidth',2)
%     hold on
%     plot([freqcomp freqcomp],[-90 90],'-k','LineWidth',2)
%     hold on
%     subplot(614);
%     title(sprintf('1-2, Station: %s',sta))
%     xlabel('Frequency (Hz)')
%     ylabel('Phase')
%     ylim([-90 90])
%     plot([pper pper],[-90 90],'-b','LineWidth',2)
%     hold on
%     plot([sper sper],[-90 90],'-r','LineWidth',2)
%     hold on
%     plot([freqcomp freqcomp],[-90 90],'-k','LineWidth',2)
%     hold on
%     subplot(615);
%     title(sprintf('1-P, Station: %s',sta))
%     xlabel('Frequency (Hz)')
%     ylabel('Phase')
%     ylim([-90 90])
%     plot([pper pper],[-90 90],'-b','LineWidth',2)
%     hold on
%     plot([sper sper],[-90 90],'-r','LineWidth',2)
%     hold on
%     plot([freqcomp freqcomp],[-90 90],'-k','LineWidth',2)
%     hold on
%     subplot(616);
%     title(sprintf('2-P, Station: %s',sta))
%     xlabel('Frequency (Hz)')
%     ylabel('Phase')
%     ylim([-90 90])
%     plot([pper pper],[-90 90],'-b','LineWidth',2)
%     hold on
%     plot([sper sper],[-90 90],'-r','LineWidth',2)
%     hold on
%     plot([freqcomp freqcomp],[-90 90],'-k','LineWidth',2)
%     hold on
%     set(gcf,'PaperPositionMode','manual');
%     set(gcf,'PaperUnits','inches');
%     set(gcf,'PaperOrientation','portrait');
%     set(gcf,'PaperPosition',[.05 .05 8 10.5]);
%     filename=sprintf('%s/Phase_%s',figoutpath,sta);
%     print(gcf,'-dpng',filename)
% 
%     figure(5)
%     subplot(611);
%     title(sprintf('Z-1, Station: %s',sta))
%     xlabel('Frequency (Hz)')
%     ylabel('Admittance')
%     hold on
%     subplot(612);
%     title(sprintf('Z-2, Station: %s',sta))
%     xlabel('Frequency (Hz)')
%     ylabel('Admittance')
%     hold on
%     subplot(613);
%     title(sprintf('P-Z, Station: %s',sta))
%     xlabel('Frequency (Hz)')
%     ylabel('Admittance')
%     hold on
%     subplot(614);
%     title(sprintf('1-2, Station: %s',sta))
%     xlabel('Frequency (Hz)')
%     ylabel('Admittance')
%     hold on
%     subplot(615);
%     title(sprintf('1-P, Station: %s',sta))
%     xlabel('Frequency (Hz)')
%     ylabel('Admittance')
%     hold on
%     subplot(616);
%     title(sprintf('2-P, Station: %s',sta))
%     xlabel('Frequency (Hz)')
%     ylabel('Admittance')
%     hold on
%     set(gcf,'PaperPositionMode','manual');
%     set(gcf,'PaperUnits','inches');
%     set(gcf,'PaperOrientation','portrait');
%     set(gcf,'PaperPosition',[.05 .05 8 10.5]);
%     filename=sprintf('%s/Admittance_%s',figoutpath,sta);
%     print(gcf,'-dpng',filename)



