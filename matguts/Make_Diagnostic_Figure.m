cmappp2plt = jet(length(indgd));
Diagnose=0;
FreqAx=[0.01 2.5];

figmetric=figure(200);
subplot(2,3,[2 5])
spectsum=[];
for ijklmn = 1:length(indgd)
currdx = indgd(ijklmn);
temptempttt = tts;
temptemptdat = all_datc(:,currdx);
dt = temptempttt(2)-temptempttt(1);
npass=2;
% loop over filterbank

for iijjkkllmm = 1:length(flos)
flo=   flos(iijjkkllmm);
fhi=fhis(iijjkkllmm);
centfreq(iijjkkllmm)=(fhi+flo)/2;
recf = filt_quick(temptemptdat,flo,fhi,dt,cp.npoles,npass);

WfStore(ijklmn,iijjkkllmm,1:length(temptempttt)) = recf./max(recf)+iijjkkllmm;
plot(temptempttt,recf./max(recf)+iijjkkllmm,'linewidth',0.5,'color',[cmappp2plt(ijklmn,:) 0.25])
hold on

sigwin = find(temptempttt > -1/centfreq(iijjkkllmm) & temptempttt < 1.5/centfreq(iijjkkllmm));
noisewin = find(temptempttt < -1/centfreq(iijjkkllmm) & temptempttt > -3/centfreq(iijjkkllmm) ); % only use time window before P wave because it's the fastest phase
snrwf(ijklmn,iijjkkllmm) = (rms(recf(sigwin))./rms(recf(noisewin)))^2;


end
Fs = 1/dt;
TW=(FreqResolution*WinSize)/2;
NTaper=2*TW-1;
taper_params = [TW NTaper];
window_params =[WinSize StepSize];

[spect,stimes,sfreqs] = multitaper_spectrogram(temptemptdat,...
Fs, FreqAx, taper_params, window_params,...
0, 'linear', 'unity', false, false);


nan_dB = nanpow2db(spect);

if ijklmn == 1
spectsum = nan_dB;
else
spectsum=[spectsum + nan_dB];
end


end

figure(200)
subplot(2,3,[1 4])
h= pcolor(stimes-cp.pretime, sfreqs, (spectsum));
h.EdgeColor = 'none';
xlim([-40 40])
cmap=crameri('imola');
colormap(cmap)
axis xy
barbar=colorbar;
ylabel(barbar,'db')
title('Stacked Spectrogram')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gca,'fontsize',18,'YScale','log')


% Now plot avg. freq-dependent snr        
figure(200)

subplot(2,3,6)
for abc = 1:length(flos)
plot([flos(abc) fhis(abc)],[mean(snrwf(:,abc))  mean(snrwf(:,abc))],'linewidth',4,'color','k') 
hold on 
errlist(abc) = std(snrwf(:,abc));
meansnrlist(abc) = mean(snrwf(:,abc));
end
xlabel('Center Frequency (Hz)')
ylabel('Signal-to-Noise Ratio')
set(gca,'fontsize',16,'fontweight','bold')
grid on; box on; grid minor;
ghj=errorbar(centfreq,meansnrlist,errlist,'linewidth',2,'color','k');
ghj.CapSize=18;
set(gca,'XScale','log')

figure(200)

subplot(2,3,[2 5])
xlabel('Time (s)')
for iijjkkllmm = 1:length(flos)
flo=   flos(iijjkkllmm);
fhi=fhis(iijjkkllmm);
centfreq(iijjkkllmm)=(fhi+flo)/2;
text(-44, iijjkkllmm+0.5,[num2str(round(flo,2)) '-' num2str(round(fhi,2)) ' Hz, ' num2str(round(1/flo,2)) '-' num2str(round(1/fhi,2)) 's'],'fontsize',17)  
set(gca,'YTickLabel',[]);
ylim([0 length(flos)+1])
end
xlim([-50 20])
set(gca,'fontsize',16,'fontweight','bold')


stack_wf = sum(all_datc(:,indgd)')./length(indgd);

for iijjkkllmm = 1:length(flos)
flo=   flos(iijjkkllmm);
fhi=fhis(iijjkkllmm);
centfreq(iijjkkllmm)=(fhi+flo)/2;
filtstack = filt_quick(stack_wf',flo,fhi,dt,cp.npoles,npass);
%filtstack=filtstack./max(filtstack);
timezoom = find(temptempttt > -1.5/centfreq(iijjkkllmm) & temptempttt < 2/centfreq(iijjkkllmm));

for asdf = 1:1:length(indgd)
   % asdf
currdx = indgd(asdf);
temptemptdat2 = all_datc(:,currdx);
filtdat = filt_quick(temptemptdat2,flo,fhi,dt,cp.npoles,npass);
%filtdat=filtdat./max(filtdat);
wintmp = tukeywin(length(timezoom),0.25);
c = xcorr(filtstack(timezoom).*wintmp,filtdat(timezoom).*wintmp,'normalized');
maxcorr = max(c);
Corrlist(asdf,iijjkkllmm) = maxcorr;

if Diagnose == 1 && iijjkkllmm == 1
    figure(99)
plot(temptempttt(timezoom),(filtstack(timezoom).*wintmp)./max(abs(filtstack(timezoom)))+asdf,'linewidth',1,'color','k')
hold on
plot(temptempttt(timezoom),(filtdat(timezoom).*wintmp)./max(abs(filtdat(timezoom)))+asdf,'linewidth',2,'color','b')

end


end
meancorrlist(iijjkkllmm) = mean(Corrlist(:,iijjkkllmm));
stdcorrlist(iijjkkllmm) = std(Corrlist(:,iijjkkllmm));


end
figmetric=figure(200);

subplot(2,3,[3])
for abc = 1:length(flos)
plot([flos(abc) fhis(abc)],[meancorrlist(abc)  meancorrlist(abc)],'linewidth',4,'color','k') 
hold on 
end
hhhhhh=errorbar(centfreq,meancorrlist,meancorrlist,'linewidth',2,'color','k');
hhhhhh.CapSize=18;
set(gca,'XScale','log')
ylim([0 1])
xlabel('Center Frequency (Hz)')
ylabel('Average Stack/Wf XCorr Coefficient')
set(gca,'fontsize',16,'fontweight','bold')
grid on; box on; grid minor;
set(figmetric,'position',[320 334 1756 980])
