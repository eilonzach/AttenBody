function [ ofile ] = plot_quick_spectra( Zraw,H1raw,H2raw,Praw,tt,Nwin,goodwin,ifsave )
% [ ofile ] = plot_quick_spectra( Zraw,H1raw,H2raw,Praw,tt,Nwin,goodwin,ifsave )
% 
% quick script to plot data or noise spectra over a number of windows:
% splits the data up into Nwin windows. Only keeps & plots the good
% windows (by default, all windows are 'good')
% 
% cribbed from script by H. Janiszewski
% adapted by Z. Eilon 03/2016

if nargin<7
    goodwin = [1:Nwin];
end
if nargin<8
    ifsave = 0;
end

%% prelims
tt = tt-min(tt);
dt = tt(2)-tt(1);
samprate = 1./dt;
npts = length(tt);
overlap = ceil((length(tt)*Nwin - length(Zraw))/(Nwin - 1));
taperf = overlap/npts;

%% compute spectra
[spect_Z,FZ,TZ] = spectrogram(Zraw,flat_hanning(tt,taperf*dt*npts),overlap,npts,samprate,'yaxis');
[spect_1,F1,T1] = spectrogram(H1raw,flat_hanning(tt,taperf*dt*npts),overlap,npts,samprate,'yaxis');
[spect_2,F2,T2] = spectrogram(H2raw,flat_hanning(tt,taperf*dt*npts),overlap,npts,samprate,'yaxis');
[spect_P,FP,TP] = spectrogram(Praw,flat_hanning(tt,taperf*dt*npts),overlap,npts,samprate,'yaxis');

spect_Z = spect_Z(:,goodwin);
spect_1 = spect_1(:,goodwin);
spect_2 = spect_2(:,goodwin);
spect_P = spect_P(:,goodwin);

TZ = TZ(goodwin);
T1 = T1(goodwin);
T2 = T2(goodwin);
TP = TP(goodwin);


%% normalise
[HZ1min,HZ1in] = min(abs(FZ-1));
czmax = max(max(log(abs(spect_Z(1:HZ1in,:)))));
czmin = min(min(log(abs(spect_Z(1:HZ1in,:)))));
c1max = max(max(log(abs(spect_1(1:HZ1in,:)))));
c1min = min(min(log(abs(spect_1(1:HZ1in,:)))));
c2max = max(max(log(abs(spect_2(1:HZ1in,:)))));
c2min = min(min(log(abs(spect_2(1:HZ1in,:)))));
cpmax = max(max(log(abs(spect_P(1:HZ1in,:)))));
cpmin = min(min(log(abs(spect_P(1:HZ1in,:)))));


%% plot
f96 = figure(96);
clf
subplot(411)
surf(TZ,FZ,log(abs(spect_Z)),'EdgeColor','none')
ylim([0 .5]);
xlim([min(TZ) max(TZ)]);
title('Spectrogram Z')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
view(2)
caxis([czmin,czmax])
colorbar
grid off
subplot(412)
surf(T1,F1,log(abs(spect_1)),'EdgeColor','none')
ylim([0 .5]);
xlim([min(T1) max(T2)]);
title('Spectrogram H1')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
view(2)
caxis([c1min,c1max])
colorbar
grid off
subplot(413)
surf(T2,F2,log(abs(spect_2)),'EdgeColor','none')
ylim([0 .5]);
xlim([min(T2) max(T2)]);
title('Spectrogram H2')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
view(2)
caxis([c2min,c2max])
colorbar
grid off
subplot(414)
surf(TP,FP,log(abs(spect_P)),'EdgeColor','none')
ylim([0 .5]);
xlim([min(TP) max(TP)]);
title('Spectrogram P')
xlabel('Time (s)')
view(2)
caxis([cpmin,cpmax])
colorbar
grid off
ylabel('Frequency (Hz)')
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);

%% save file
if ifsave
ofile = input('Give file name for figure with spectra: ','s');
print(gcf,'-dpng',ofile)
else
    ofile = [];
end


end

