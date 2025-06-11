%
Tmids = (Tlos+This)./2;
flos = 1./This;
fhis = 1./Tlos;
cfreq = (fhis + flos)/2;

% for plotting
f2plt = sort([5,4,2,1,1./[2 3 4 6 8 10 12 15 20 25 30 40 50]]');

cmappp2plt = brewermap(length(indgd),'PuBuGn');
Diagnose=0;
FreqAx=[0.01 2.5];

% setup
figmetric=figure(200); clf, set(figmetric,'pos',[-2014 241 1298 709])
axfm1 = axes('pos',[0.3808 0.100 0.25 0.8150]); hold on
axfm2 = axes('pos',[0.0800 0.1100 0.25 0.8150]); hold on
axfm3 = axes('pos',[0.6926 0.1100 0.26 0.3412]); hold on
% axfm4 = subplot(2,3,3);

spectsum=[];

%% loop over all good traces and make spectrograms
WfStore = nan(length(indgd),length(flos),length(tts));
for ig = 1:length(indgd)
is = indgd(ig);
tempttt = tts;
temptdat = all_datc(:,is);
dt = tempttt(2)-tempttt(1);

npass=2;

% for this phase, loop over filterbank
for ifilt = 1:length(flos)

    recf = filt_quick(temptdat,flos(ifilt),fhis(ifilt),dt,cp.npoles,npass);
    
    WfStore(ig,ifilt,1:length(tempttt)) = recf./max(recf)+ifilt;
    
    plot(axfm1,tempttt,recf./max(recf)+ifilt,'linewidth',0.7,'color',[cmappp2plt(ig,:) 0.65])
    
    sigwin = find(tempttt > -1/cfreq(ifilt) & tempttt < 1.5/cfreq(ifilt));
    noisewin = find(tempttt < -1/cfreq(ifilt) & tempttt > -3/cfreq(ifilt) ); % only use time window before P wave because it's the fastest phase
    snrwf(ig,ifilt) = (rms(recf(sigwin))./rms(recf(noisewin)))^2;
end

Fs = 1/dt;
TW=(FreqResolution*WinSize)/2;
NTaper=2*TW-1;
taper_params = [TW NTaper];
window_params =[WinSize StepSize];

% make individual spectrum
[spect,stimes,sfreqs] = multitaper_spectrogram(temptdat,...
                            Fs, FreqAx, taper_params, window_params,...
                            0, 'linear', 'unity', false, false);

spec_dB_bad2nan = nanpow2db(spect);

% stack individual spectrum into total spectrum
    if ig == 1 % initialise with first spcetrum
        spectsum = spec_dB_bad2nan;
    else % add later spectra on
        spectsum = [spectsum + spec_dB_bad2nan];
    end


end

% pretty the traces plot
xlabel(axfm1,'Time (s)')
for ifilt = 1:length(flos)
    text(axfm1,-44, ifilt+0.5,[num2str(round(flos(ifilt),2)) '-' num2str(round(fhis(ifilt),2)) ' Hz, ' num2str(round(1/fhis(ifilt),2)) '-' num2str(round(1/flos(ifilt),2)) 's'],'fontsize',17)  
end
set(axfm1,'fontsize',16,'fontweight','bold','YTickLabel',[],...
    'xlim',[-50 20],'ylim',[0 length(flos)+1])


%% plot the stacked spectrogram
h_spec=pcolor(axfm2,stimes-cp.pretime, sfreqs, (spectsum));
h_spec.EdgeColor = 'none';

addpath('~/Dropbox/MATLAB/lib')
% make pretty
colormap(axfm2,brewermap(40,'-Spectral'));
axis xy
% barbar=colorbar;
% ylabel(barbar,'db')

title(axfm2,'Stacked Spectrogram')
xlabel(axfm2,'Time (s)')
ylabel(axfm2,'Frequency (Hz)')
set(axfm2,'fontsize',18,'YScale','log','xlim',[-40 40])

% add y-ticks for frequencies in use
set(axfm2,'ytick',f2plt,'yticklabel',num2str(1./f2plt))

%%  Plot avg. freq-dependent snr        
for abc = 1:length(flos)
    plot(axfm3,[flos(abc) fhis(abc)],[mean(snrwf(:,abc))  mean(snrwf(:,abc))],'linewidth',4,'color','k') 
    hold on 
%     errlist(abc) = std(snrwf(:,abc));
%     meansnrlist(abc) = mean(snrwf(:,abc));
end
xlabel(axfm3,'Center Frequency (Hz)')
ylabel(axfm3,'Signal-to-Noise Ratio')
set(axfm3,'fontsize',16,'fontweight','bold','box','on');
grid on, grid minor;
% ghj=errorbar(cfreq,meansnrlist,errlist,'linewidth',2,'color','k');
% ghj.CapSize=18;
set(axfm3,'XScale','log')
set(axfm3,'xtick',f2plt,'xticklabel',num2str(1./f2plt))


return

%% stack_wf = sum(all_datc(:,indgd)')./length(indgd);

for ifilt = 1:length(flos)
flo=   flos(ifilt);
fhi=fhis(ifilt);
cfreq(ifilt)=(fhi+flo)/2;
filtstack = filt_quick(stack_wf',flo,fhi,dt,cp.npoles,npass);
%filtstack=filtstack./max(filtstack);
timezoom = find(tempttt > -1.5/cfreq(ifilt) & tempttt < 2/cfreq(ifilt));

for asdf = 1:1:length(indgd)
   % asdf
is = indgd(asdf);
temptemptdat2 = all_datc(:,is);
filtdat = filt_quick(temptemptdat2,flo,fhi,dt,cp.npoles,npass);
%filtdat=filtdat./max(filtdat);
wintmp = tukeywin(length(timezoom),0.25);
c = xcorr(filtstack(timezoom).*wintmp,filtdat(timezoom).*wintmp,'normalized');
maxcorr = max(c);
Corrlist(asdf,ifilt) = maxcorr;

if Diagnose == 1 && ifilt == 1
    figure(99)
plot(tempttt(timezoom),(filtstack(timezoom).*wintmp)./max(abs(filtstack(timezoom)))+asdf,'linewidth',1,'color','k')
hold on
plot(tempttt(timezoom),(filtdat(timezoom).*wintmp)./max(abs(filtdat(timezoom)))+asdf,'linewidth',2,'color','b')

end


end
meancorrlist(ifilt) = mean(Corrlist(:,ifilt));
stdcorrlist(ifilt) = std(Corrlist(:,ifilt));


end
figmetric=figure(200);

subplot(2,3,[3])
for abc = 1:length(flos)
plot([flos(abc) fhis(abc)],[meancorrlist(abc)  meancorrlist(abc)],'linewidth',4,'color','k') 
hold on 
end
hhhhhh=errorbar(cfreq,meancorrlist,meancorrlist,'linewidth',2,'color','k');
hhhhhh.CapSize=18;
set(gca,'XScale','log')
ylim([0 1])
xlabel('Center Frequency (Hz)')
ylabel('Average Stack/Wf XCorr Coefficient')
set(gca,'fontsize',16,'fontweight','bold')
grid on; box on; grid minor;
