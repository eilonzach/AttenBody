function plot_1evt1sta_spectra(spectra_file,ifsave,odir)
% code run after OBS_noise_cohere that reviews and plots spectral
% properties for a single station
% 
% INPUTS
%   spectra_file = .mat file with all spectra details for one station for
%                   one event, e.g. % spectra_file = '/Volumes/DATA/CASCADIA/DATA/215_201308301625/J32C_spectra.mat';
close all
if nargin < 2
    ifsave = 0;
end

if nargin < 3
    odir = './';
end

if ~strcmp(odir(end),'/'), odir = [odir,'/']; end

if ifsave
    if exist(odir,'dir')~=7, mkdir(odir); end
end

load(spectra_file);

evid = strtok(evdir,'/');

grey = [0.5 0.5 0.5];


% Plotting Power
figure(3), clf, set(gcf,'pos',[30 30 800 600])
subplot(411);
loglog(f,smooth(czz_stack,40),'-','LineWidth',.5,'Color',grey); ylabel('ZZ')
title('Power','Fontsize',16,'Interpreter','latex')
subplot(412);
loglog(f,smooth(c11_stack,40),'-','LineWidth',.5,'Color',grey); ylabel('11')
subplot(413);
loglog(f,smooth(c22_stack,40),'-','LineWidth',.5,'Color',grey); ylabel('22')
subplot(414);
loglog(f,smooth(cpp_stack,40),'-','LineWidth',.5,'Color',grey); ylabel('PP')
if ifsave
    filename = [odir,sta,'_',evid,'_power'];
    print(gcf,'-dpng',filename)
end

% Plotting Coherence
figure(4), clf, set(gcf,'pos',[40 40 800 600])
subplot(611);
semilogx(f,smooth(coh1z_stack,40),'-','LineWidth',.5,'Color',grey); ylabel('1Z')
title('Coherence','Fontsize',16,'Interpreter','latex')
subplot(612);
semilogx(f,smooth(coh2z_stack,40),'-','LineWidth',.5,'Color',grey); ylabel('2Z')
subplot(613);
semilogx(f,smooth(cohpz_stack,40),'-','LineWidth',.5,'Color',grey); ylabel('PZ')
subplot(614);
semilogx(f,smooth(coh12_stack,40),'-','LineWidth',.5,'Color',grey); ylabel('12')
subplot(615);
semilogx(f,smooth(coh1p_stack,40),'-','LineWidth',.5,'Color',grey); ylabel('1P')
subplot(616);
semilogx(f,smooth(coh2p_stack,40),'-','LineWidth',.5,'Color',grey); ylabel('2P')
if ifsave
    filename = [odir,sta,'_',evid,'_coherence'];
    print(gcf,'-dpng',filename)
end

% Plotting Phase
figure(6), clf, set(gcf,'pos',[60 60 800 600])
subplot(611);
semilogx(f,ph1z_stack,'-','LineWidth',.5,'Color',grey); ylabel('1Z')
title('Phase','Fontsize',16,'Interpreter','latex')
subplot(612);
semilogx(f,ph2z_stack,'-','LineWidth',.5,'Color',grey); ylabel('2Z')
subplot(613);
semilogx(f,phpz_stack,'-','LineWidth',.5,'Color',grey); ylabel('PZ')
subplot(614);
semilogx(f,ph12_stack,'-','LineWidth',.5,'Color',grey); ylabel('12')
subplot(615);
semilogx(f,ph1p_stack,'-','LineWidth',.5,'Color',grey); ylabel('1P')
subplot(616);
semilogx(f,ph2p_stack,'-','LineWidth',.5,'Color',grey); ylabel('2P')
if ifsave
    filename = [odir,sta,'_',evid,'_phase'];
    print(gcf,'-dpng',filename)
end


% Plotting Admittance
figure(5), clf, set(gcf,'pos',[50 50 800 600])
subplot(611);
loglog(f,smooth(ad1z_stack,40),'-','LineWidth',.5,'Color',grey); ylabel('1Z')
title('Admittance','Fontsize',16,'Interpreter','latex')
subplot(612);
loglog(f,smooth(ad2z_stack,40),'-','LineWidth',.5,'Color',grey); ylabel('2Z')
subplot(613);
loglog(f,smooth(adpz_stack,40),'-','LineWidth',.5,'Color',grey); ylabel('PZ')
subplot(614);
loglog(f,smooth(ad12_stack,40),'-','LineWidth',.5,'Color',grey); ylabel('12')
subplot(615);
loglog(f,smooth(ad1p_stack,40),'-','LineWidth',.5,'Color',grey); ylabel('1P')
subplot(616);
loglog(f,smooth(ad2p_stack,40),'-','LineWidth',.5,'Color',grey); ylabel('2P')
if ifsave
    filename = [odir,sta,'_',evid,'_admittance'];
    print(gcf,'-dpng',filename)
end

