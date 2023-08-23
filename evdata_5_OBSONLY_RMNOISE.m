% Script to compute spectral characteristics of noise before each event.
% This script goes through data event by event, then station by station.
% For OBS stations, it pulls out the 12 hours of noise prior to each event
% and computes spectra on a set of Nwin windows of noise, saving the
% outputs into structures within the data directories
clear all
close all


% PICK NOISE CORRECTION METHOD
% different methods to compute transfer function - from Helen
% 1 = TILT AND COMPLIANCE, CRAWFORD AND WEBB, 2000, spahr's code
% 2 = COMPLIANCE ONLY, same nomenclature but only removing pressure component
% 3 = COMPLIANCE AND THEN TILT, individual component method
% (4) = TILT AND COMPLIANCE - BELL ET AL., 2014 reorient in maximum coherence direction
method_trans_func = 2; 

% CHOOSE WHICH DATA TO CORRECT
cor_opt = 1; % option for which data to correct: 0) all, 1) data, 2) noise
 % time windows for data excerpt - from event time
event_win = [-100 1700]; % can be any length, but refer to evdata_1_DOWNLOAD
noise_win = [-20000 -14000]; % should be T seconds (6000) in length, refer to evdata_4_NOISESPECTRA

% OPTIONS TO SMOOTH TRANSFER FUNCTIONS
smooth_trans_func = false; % option to smooth transfer functions
npt_mav_stf = 10; %moving average window```

% OPTION TO FILTER TRANSFER FUNCTIONS - definitely do this!!
filt_trans_func = true; % option to low-pass filter transfer functions 

% OTHER OPTIONS
overwrite = true;
ifplot = 0; % option to do the plots
ifsavefigs = 0; % option to save the figures
ifmanualQC = 0; % option to manually choose if to keep corrected data

% RESAMPLE DATA - ## CURRENTLY NOT CODED ##
resamprate = 5; % new sample rate to downsamp to

% SPECIFY SOME BASIC PARMS
pper = 1/16;
sper = 1/8;

fmin=0.005;fmax=0.02;
coh_min=0.8;
chorz_max=1;

plot_filter = [1 50]; %[lowpass highpass] in seconds (i.e. [min_period max_period])

periods = [1 10 20 50 75 100 125 150 200];


%% paths
cd('/Users/zeilon/Documents/MATLAB/CASC_atten/')
addpath('/Users/zeilon/Documents/MATLAB/helen_dataprocessing')
addpath('matguts')
% path to sensor spectra data
specdir = '~/Work/CASCADIA/CAdb/spectra/';
% antelope db details
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% path to top level of directory tree for data
datadir = '/Volumes/DATA_mini2/CASCADIA/DATA/'; % needs final slash


%% get to work
[ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam ); % load events data
[ nstas,stas,slats,slons,selevs,ondate,offdate,staname,statype ] = db_stadata( dbdir,dbnam ); %load stations data

for ie = 1:269 % 1:norids % loop on orids
    fprintf('\n Orid %.0f %s \n\n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
    evdir = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
    evday = epoch2str(evtimes(ie),'%Y%j');
    datinfofile = [datadir,evdir,'_datinfo'];
    evstr = regexprep(strtok(evdir,'/'),'_','-');
           
    if any((evtimes(ie)-evtimes)>0 & (evtimes(ie)-evtimes) < 20*60*60)
        fprintf('Another event within prev 20 hrs... skipping\n')
        continue
    end

    % check files exist
    if exist([datinfofile,'.mat'],'file')~=2, fprintf('No data at all for this event\n');continue, end
    load(datinfofile)
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end

    for is = 1:length(datinfo) % loop on stas
        sta = datinfo(is).sta; % sta name
        
        %% load dataf, check for previous activity
        if strcmp(statype(strcmp(stas,sta)),'OBS') ~= 1, continue, end % skip if not an OBS
        fprintf('Station %-3.0f, %-5s...',is,sta)
        
        specfile = [datadir,evdir,sta,'_spectra.mat'];
        datfile = [datadir,evdir,sta,'.mat'];
        
        % check files exist, otherwise skip
        if exist(specfile,'file')~=2, fprintf(' no spectra file.\n'), continue, end
        if exist(datfile,'file')~=2, fprintf(' no data file.\n'),  continue, end
        % conditions before loading - don't unecessarily load
        if datinfo(is).rmtilt && datinfo(is).rmtilt && ~overwrite
            fprintf(' done already.\n'), continue
        end

        % OKAY - load data
        load(datfile)
        load(specfile)
        
        % second pass on conditions
        if data.rmtilt && data.rmcomp && ~overwrite, fprintf(' done already.\n'), continue, end % skip if already removed tilt + comp
        
        %% calc some basics
        T = (win_pt_end(1)-win_pt_start(1)+1)*dt;
        Nwin = length(win_pt_start);
        samprate = 1./dt;
        freqcomp = sqrt(9.8/(-2*pi*data.station.selev));

%% ==================================================================== %%       
        %% COMPUTE TRANSFER FUNCTION
%% ==================================================================== %%   
        fprintf(' calc transfunc...')
        
        switch method_trans_func
            case 1
                [czz_12p,czz_12,czz_1,cohpz_12,l1z,l12,l1p,l2p_1,l2z_1,lpz_12,coh2z_1,cohpz_1,coh2p_1]...
                    = multicoher(c11_stack,c22_stack,cpp_stack,czz_stack,c1z_stack,c2z_stack,cpz_stack,c12_stack,c1p_stack,c2p_stack,f);
                transfer_fun = struct('sta',sta,'method',1,'f',f,'freqcomp',freqcomp,'NFFT',NFFT,...
                                    'czz_12p',czz_12p,'czz_12',czz_12,'czz_1',czz_1,...
                                    'l1z',l1z,'l12',l12,'l1p',l1p,'l2p_1',l2p_1,'l2z_1',l2z_1,'lpz_12',lpz_12,...
                                    'cohpz_12',cohpz_12,'coh2z_1',coh2z_1,'cohpz_1',cohpz_1,'coh2p_1',coh2p_1);
            case 2
                [czz_p,lpz,lp1,lp2,coh1z_p,coh2z_p,coh12_p,c11_p,c22_p,c1z_p,c2z_p,c12_p]...
                    = compcohere(c11_stack,c22_stack,cpp_stack,czz_stack,c1z_stack,c2z_stack,cpz_stack,c12_stack,c1p_stack,c2p_stack,f);
                transfer_fun = struct('sta',sta,'method',2,'f',f,'freqcomp',freqcomp,'NFFT',NFFT,...
                                    'czz_p',czz_p,'lpz',lpz,'lp1',lp1,'lp2',lp2,...
                                    'coh1z_p',coh1z_p,'coh2z_p',coh2z_p,'coh12_p',coh12_p,...
                                    'c11_p',c11_p,'c22_p',c22_p,'c1z_p',c1z_p,'c2z_p',c2z_p,'c12_p',c12_p);
            case 3
                [czz_p21,czz_p2,czz_p,coh1z_p2,lpz,lp2,lp1,l21_p,l2z_p,l1z_p2,coh2z_p,coh1z_p,coh21_p]...
                    = multicoher(cpp_stack,c22_stack,c11_stack,czz_stack,cpz_stack,c2z_stack,c1z_stack,c2p_stack,c1p_stack,c12_stack,f);
                transfer_fun = struct('sta',sta,'method',3,'f',f,'freqcomp',freqcomp,'NFFT',NFFT,...
                                    'czz_p21',czz_p21,'czz_p2',czz_p2,'czz_p',czz_p,...
                                    'lpz',lpz,'lp1',lp1,'lp2',lp2,'l21_p',l21_p,'l2z_p',l2z_p,'l1z_p2',l1z_p2,...
                                    'coh1z_p2',coh1z_p2,'coh2z_p',coh2z_p,'coh1z_p',coh1z_p,'coh21_p',coh21_p);
            case 4
                fprintf('Sorry, this method not yet ready for primetime!\n'), return
        end
        
        %% SMOOTH THE TRANSFER FUNCTIONS
        if smooth_trans_func
            switch method_trans_func
                case 1
                    % Filter transfer function to avoid microseism band
                    l1z=smooth(l1z,npt_mav_stf);
                    l12=smooth(l12,npt_mav_stf);
                    l1p=smooth(l1p,npt_mav_stf);
                    l2p_1=smooth(l2p_1,npt_mav_stf);
                    l2z_1=smooth(l2z_1,npt_mav_stf);
                    lpz_12=smooth(lpz_12,npt_mav_stf);
                case 2
                    % Filter transfer function to avoid microseism band
                    lpz=smooth(lpz,npt_mav_stf);
                case 3
                    % Filter transfer function to avoid microseism band
                    lpz=smooth(lpz,npt_mav_stf);
                    lp2=smooth(lp2,npt_mav_stf);
                    lp1=smooth(lp1,npt_mav_stf);
                    l21_p=smooth(l21_p,npt_mav_stf);
                    l2z_p=smooth(l2z_p,npt_mav_stf);
                    l1z_p2=smooth(l1z_p2,npt_mav_stf);
            end
        end

        %% FILTER THE TRANSFER FUNCTIONS TO AVOID MICROSEISM BAND
        if filt_trans_func
            if pper < freqcomp
                lp=pper; %low pass filter transfer function at primary microseism
            else 
                lp=freqcomp;
            end        

            switch method_trans_func
                case 1
                    % Filter transfer function to avoid microseism band
                    l1z=gauss2_freq(l1z,f,lp,NFFT);
                    l12=gauss2_freq(l12,f,lp,NFFT);
                    l1p=gauss2_freq(l1p,f,lp,NFFT);
                    l2p_1=gauss2_freq(l2p_1,f,lp,NFFT);
                    l2z_1=gauss2_freq(l2z_1,f,lp,NFFT);
                    lpz_12=gauss2_freq(lpz_12,f,lp,NFFT);
                case 2
                    % Filter transfer function to avoid microseism band
                    lpz=gauss2_freq(lpz,f,lp,NFFT);
                case 3
                    % Filter transfer function to avoid microseism band
                    lpz=gauss2_freq(lpz,f,lp,NFFT);
                    lp2=gauss2_freq(lp2,f,lp,NFFT);
                    lp1=gauss2_freq(lp1,f,lp,NFFT);
                    l21_p=gauss2_freq(l21_p,f,lp,NFFT);
                    l2z_p=gauss2_freq(l2z_p,f,lp,NFFT);
                    l1z_p2=gauss2_freq(l1z_p2,f,lp,NFFT);
            end
        end

        %% SAVE smoothed/filtered transfer function back into transfer_fun structure
        transfer_fun.smoothTF = smooth_trans_func;
        transfer_fun.filtTF = filt_trans_func;
        switch method_trans_func
            case 1
                transfer_fun.l1z_use=l1z;
                transfer_fun.l12_use=l12;
                transfer_fun.l1p_use=l1p;
                transfer_fun.l2p_1_use=l2p_1;
                transfer_fun.l2z_1_use=l2z_1;
                transfer_fun.lpz_12_use=lpz_12;
            case 2
                % Filter transfer function to avoid microseism band
                transfer_fun.lpz_use=lpz;
            case 3
                % Filter transfer function to avoid microseism band
                transfer_fun.lpz_use=lpz;
                transfer_fun.lp2_use=lp2;
                transfer_fun.lp1_use=lp1;
                transfer_fun.l21_p_use=l21_p;
                transfer_fun.l2z_p_use=l2z_p;
                transfer_fun.l1z_p2_use=l1z_p2;
        end
        
        
        %% PLOT transfer functions
        if ifplot
            plot_trans_func(method_trans_func,transfer_fun) %< option in here to save, can switch on
            set(gcf,'pos',get(gcf,'pos')+[0 500 0 0])
        end
        
        %% metric for size of transfer function - i.e. is it worth removing tilt?!
        % consider looping back and doing only compliance removal if tilt
        % is little
        
%% ==================================================================== %%       
        %% CORRECT THE SEISMOGRAMS
%% ==================================================================== %%       
        fprintf(' correct Z...')

        %% GET AND PROCESS EARTHQUAKE DATA
        
        % sort out channels        
        chans = data.chans;        
        ich = find(strcmp(chans.component,'H'));
        icn = find(strcmp(chans.component,'N'));
        ice = find(strcmp(chans.component,'E')); 
        icz = find(strcmp(chans.component,'Z'));
        if isempty(icn), icn = find(strcmp(chans.component,'1')); end
        if isempty(ice), ice = find(strcmp(chans.component,'2')); end
        ics = [ich icn ice icz];
        if length(ics)~=4  % skip if don't have all 4 chans
            fprintf('only %.0f chans... skipping\n',length(ics));  
            continue, 
        end 
        
        % window to just the data we're going to correct
        if cor_opt==0
            j1 = 1;
            j2 = length(data.tt);
        elseif cor_opt==1 
            j1 = mindex(data.tt,evtimes(ie)+event_win(1));
            j2 = j1 + ( diff(event_win)*data.samprate - 1 );
        elseif cor_opt == 2
            j1 = mindex(data.tt,evtimes(ie)+noise_win(1));
            j2 = j1 + ( diff(noise_win)*data.samprate - 1 );        
        end
        
        % grab data
        Zdat  = data.dat(j1:j2,icz);
        H1dat = data.dat(j1:j2,icn);
        H2dat = data.dat(j1:j2,ice);
        Pdat  = data.dat(j1:j2,ich);
        tt    = data.tt(j1:j2)-evtimes(ie);
        
        % NEED TO USE THE RAW PRESSURE IF THERE IS A LAMONT APG
        if strcmp(which_OBS(sta),'LDEO')
        	ichr = find(strcmp(raw.chans.component,'H'));
            Pdat = data.raw.dat(j1:j2,ichr);
            lo_corner = 0.005;  % in Hz
            npoles=5;
            lo_w=2*pi*lo_corner;

            N = length(Pdat);
            delta = 1./data.samprate;
            Tr = N*delta;

            if mod(N,2)
                faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/Tr);
            else
                faxis = [0:N/2,-N/2+1:-1]*(1/Tr);
            end
            w = faxis.*2*pi;

            hpfiltfrq=( ((w./lo_w).^(2*npoles))./(1+(w./lo_w).^(2*npoles)) );
            norm_trans=hpfiltfrq;    % this is normalization transfer function
            norm_trans(isnan(norm_trans)) = 0;

            fftdata = fft(Pdat);
            fftdata = fftdata(:).*norm_trans(:);
            Pdat = real(ifft(fftdata));
        end
        
        % downsamp here!
        if resamprate~=data.samprate
            fprintf('NEED TO CODE UP DATA RESAMP!\n')
            % put resampling in here (just compute new tt vector and
            % interpolate the data onto it)
        end
        
        %% PLOT EVENT DATA to check all components are working during event time
        fNq = 1/2/dt;
        T1 = plot_filter(1); T2 = plot_filter(2);
        [b,a]=butter(2,[1/fNq/T2,1/fNq/T1]);

        Z_filt  = filtfilt(b,a,Zdat);
        H1_filt  = filtfilt(b,a,H1dat);
        H2_filt  = filtfilt(b,a,H2dat);
        P_filt  = filtfilt(b,a,Pdat);

        if ifplot % plot original seismograms
            figure(101), clf; set(gcf,'pos',[-1200 340 1100 800])
            subplot(411)
            title(sprintf('%s  %s',sta,evstr)), hold on
            plot(tt,Z_filt,'-k'); ylabel('Z')
            subplot(412)
            plot(tt,H1_filt,'-k'); ylabel('H1')
            subplot(413)
            plot(tt,H2_filt,'-k'); ylabel('H2')
            subplot(414)
            plot(tt,P_filt,'-k'); ylabel('P')

            set(gcf,'PaperPositionMode','manual');
            set(gcf,'PaperUnits','inches');
            set(gcf,'PaperOrientation','portrait');
            set(gcf,'PaperPosition',[.05 .05 10.5 8]);  
        end
        

        %% COMPUTE DATA SPECTRA
        % calc. various parameters
        npts = length(Zdat);
        Ppower = nextpow2(npts);
        NFFT = 2^Ppower;
        % NFFT = npts;
        samprate = 1/dt; %common named Fs
        fdat = samprate/2*linspace(0,1,NFFT/2+1);
        npad0 = floor((NFFT-npts)/2);

        amp_Z  = Zdat;
        amp_Z  = amp_Z.*flat_hanning(tt,0.1*dt*npts);
        amp_Z  = detrend(amp_Z,0);
        amp_Z  = padarray(amp_Z,[npad0 0],'both');
        spectrum = fft(amp_Z,NFFT).*dt;
        spec_Z = spectrum(1:NFFT/2+1);

        amp_H1  = H1dat;
        amp_H1  = amp_H1.*flat_hanning(tt,0.1*dt*npts);
        amp_H1  = detrend(amp_H1,0);
        amp_H1  = padarray(amp_H1,[npad0 0],'both');
        spectrum = fft(amp_H1,NFFT).*dt;
        spec_H1 = spectrum(1:NFFT/2+1);

        amp_H2  = H2dat;
        amp_H2  = amp_H2.*flat_hanning(tt,0.1*dt*npts);
        amp_H2  = detrend(amp_H2,0);
        amp_H2  = padarray(amp_H2,[npad0 0],'both');
        spectrum = fft(amp_H2,NFFT).*dt;
        spec_H2 = spectrum(1:NFFT/2+1);

        amp_P  = Pdat;
        amp_P  = amp_P.*flat_hanning(tt,0.1*dt*npts);
        amp_P  = detrend(amp_P,0);
        amp_P  = padarray(amp_P,[npad0 0],'both');
        spectrum = fft(amp_P,NFFT).*dt;
        spec_P = spectrum(1:NFFT/2+1);

%%  ------------------------- REMOVE THE NOISE!! ------------------------
        

        %% do noise removal
        switch method_trans_func
            case 1 % =======  TILT AND COMPLIANCE METHOD 1  ======= 
                % interpolate transfer functions to the frequencies for the data
                l1z_use = interp1(transfer_fun.f,l1z,fdat);
                l12_use = interp1(transfer_fun.f,l12,fdat);
                l2z_1_use = interp1(transfer_fun.f,l2z_1,fdat);
                l1p_use = interp1(transfer_fun.f,l1p,fdat);
                l2p_1_use = interp1(transfer_fun.f,l2p_1,fdat);
                lpz_12_use = interp1(transfer_fun.f,lpz_12,fdat);

                specZ_1 = spec_Z-((l1z_use)'.*spec_H1); % remove 1 from Z ==> Zcor1
                spec2_1 = spec_H2-(l12_use'.*spec_H1); % remove 1 from 2 ==> 2cor1
                specZ_12 = specZ_1-(l2z_1_use'.*spec2_1); % remove 2cor1 from Zcor1 ==> Zcor12
                specP_1 = spec_P-(l1p_use'.*spec_H1); % remove 1 from P ==> Pcor1
                specP_12 = specP_1-(l2p_1_use'.*spec2_1); % remove 2cor1 from Pcor1 ==> Pcor12
                specZ_12P = specZ_12-(lpz_12_use'.*specP_12); % remove Pcor12 from Zcor12 => Zcor12P

                % power spectra
                poweror = abs(spec_Z).^2*2/(NFFT*dt);
                powerco = abs(specZ_12).^2*2/(NFFT*dt);
                powerc2 = abs(specZ_12P).^2*2/(NFFT*dt);
                
                % transform back to time domain
                ampZ_1 = real(ifft(2*specZ_1,NFFT)./dt);
                ampZ_12 = real(ifft(2*specZ_12,NFFT)./dt);
                ampZ_12P = real(ifft(2*specZ_12P,NFFT)./dt);
                amp_Z_1_nopad = ampZ_1(npad0+1:npad0+npts); %Raw Z data -1
                amp_Z_12_nopad = ampZ_12(npad0+1:npad0+npts); %Raw Z data -12
                amp_Z_12P_nopad = ampZ_12P(npad0+1:npad0+npts); %Raw Z data -12P
                Z_filt_1 = filtfilt(b,a,amp_Z_1_nopad);
                Z_filt_12 = filtfilt(b,a,amp_Z_12_nopad);
                Z_filt_12P = filtfilt(b,a,amp_Z_12P_nopad);

                % final corrected data
                Zdat_cor = Z_filt_12P;
                
                %% Some plots
                if ifplot
                % -------- plot progressive corrections to the Z comp ------
                figure(102), clf, set(gcf,'pos',[20 20 1100 700])
                subplot(411)
                title(sprintf('%s  %s',sta,evstr)), hold on
                plot(tt,Z_filt,'-k'); ylabel('Z')
                subplot(412)
                plot(tt,Z_filt_1,'-k'); ylabel('Z-1')
                subplot(413)
                plot(tt,Z_filt_12,'-k'); ylabel('Z-12')
                subplot(414)
                plot(tt,Z_filt_12P,'-k'); ylabel('Z-12P')
                
                set(gcf,'PaperPositionMode','manual');
                set(gcf,'PaperUnits','inches');
                set(gcf,'PaperOrientation','portrait');
                set(gcf,'PaperPosition',[.05 .05 10.5 8]);
                
                % -------- plot improvements in narrow filter bands ------
                figure(103), clf, set(gcf,'pos',[1050 320 1000 1000])
                for j = 1:length(periods)
                    T1p = periods(j)-0.1*periods(j); T2p= periods(j)+0.1*periods(j);
                    [bp,ap]=butter(2,[1/fNq/T2p,1/fNq/T1p]);
                    Z_filtp  = filtfilt(bp,ap,Zdat);
                    Z_12_filtp = filtfilt(bp,ap,amp_Z_12_nopad);
                    Z_12P_filtp = filtfilt(bp,ap,amp_Z_12P_nopad);

                    subplot(length(periods),3,3*j-2)
                    plot(tt,(Z_filtp),'-k');
                    title(sprintf('orid %s   Z   %d s',strtok(evstr,'-'), periods(j)));

                    subplot(length(periods),3,(3*j)-1)
                    plot(tt,(Z_12_filtp),'-k');
                    title(sprintf('orid %s   Z-12   %d s',strtok(evstr,'-'), periods(j)));

                    subplot(length(periods),3,(3*j))
                    plot(tt,(Z_12P_filtp),'-k');
                    title(sprintf('orid %s   Z-12P   %d s',strtok(evstr,'-'), periods(j)));
                end
                
                set(gcf,'PaperPositionMode','manual');
                set(gcf,'PaperUnits','inches');
                set(gcf,'PaperOrientation','portrait');
                set(gcf,'PaperPosition',[.05 .05 11 12]);
                               
                % -------- plot power spectra through correction ------                
                figure(104), clf, set(gcf,'pos',[2000 850 550 480]), hold on
                loglog(fdat,smooth(poweror,40),'-k','LineWidth',.5)
                loglog(fdat,smooth(powerco,40),'-r','LineWidth',.5)
                loglog(fdat,smooth(powerc2,40),'-b','LineWidth',.5)
                legend('Original','RM-tilt','RM-comp')

                title(sprintf('Z-component, Event: %s',evstr));
                xlabel('Frequency (Hz)')
                ylabel('Power Spectrum')
                set(gca,'xscale','log','yscale','log')
                
                set(gcf,'PaperPositionMode','manual');
                set(gcf,'PaperUnits','inches');
                set(gcf,'PaperOrientation','portrait');
                set(gcf,'PaperPosition',[.05 .05 7 5]);

                end % ifplot
            
            case 2 % =======  COMPLIANCE CORRECTION ONLY  ======= 
                % interpolate transfer function to the frequencies for the data
                lpz_use = interp1(transfer_fun.f,lpz,fdat);

                spectrumZ_P = spec_Z-((lpz_use').*spec_P); % remove P from Z ==> ZcorP
        
                % power spectra
                poweror = abs(spec_Z).^2*2/(NFFT*dt);
                powerco = abs(spectrumZ_P).^2*2/(NFFT*dt);
                
                % transform back to time domain
                ampZ_P = real(ifft(2*spectrumZ_P,NFFT)./dt);
                ampZ_P_nopad = ampZ_P(npad0+1:npad0+npts); %Raw Z data -P
                Z_filt_P = filtfilt(b,a,ampZ_P_nopad);

                % final corrected data
                Zdat_cor = Z_filt_P;

                %% Some plots
                if ifplot
                % -------- plot progressive corrections to the Z comp ------
                figure(102), clf, set(gcf,'pos',[20 20 1100 500])
                subplot(211)
                title(sprintf('%s  %s',sta,evstr)), hold on
                plot(tt,Z_filt,'-k'); ylabel('Z')
                subplot(212)
                plot(tt,Z_filt_P,'-k'); ylabel('Z-P')
                
                set(gcf,'PaperPositionMode','manual');
                set(gcf,'PaperUnits','inches');
                set(gcf,'PaperOrientation','portrait');
                set(gcf,'PaperPosition',[.05 .05 10.5 8]);
                
                % -------- plot improvements in narrow filter bands ------
                figure(103), clf, set(gcf,'pos',[1050 150 800 1000])
                for j = 1:length(periods)
                    T1p = periods(j)-0.1*periods(j); T2p= periods(j)+0.1*periods(j);
                    [bp,ap]=butter(2,[1/fNq/T2p,1/fNq/T1p]);
                    Z_filtp  = filtfilt(bp,ap,Zdat);
                    Z_P_filtp = filtfilt(bp,ap,ampZ_P_nopad);

                    subplot(length(periods),2,2*j-1)
                    plot(tt,(Z_filtp),'-k');
                    title(sprintf('orid %s   Z   %d s',strtok(evstr,'-'), periods(j)));

                    subplot(length(periods),2,2*j)
                    plot(tt,(Z_P_filtp),'-k');
                    title(sprintf('orid %s   Z-P   %d s',strtok(evstr,'-'), periods(j)));
                end

                set(gcf,'PaperPositionMode','manual');
                set(gcf,'PaperUnits','inches');
                set(gcf,'PaperOrientation','portrait');
                set(gcf,'PaperPosition',[.05 .05 11 12]);

                                
                % -------- plot power spectra through correction ------                
                figure(104), clf, hold on
                loglog(fdat,smooth(poweror,40),'-k','LineWidth',.5)
                loglog(fdat,smooth(powerco,40),'-r','LineWidth',.5)
                legend('Original','RM-comp')

                title(sprintf('Z-component, Event: %s',evstr));
                xlabel('Frequency (Hz)')
                ylabel('Power Spectrum')
                set(gca,'xscale','log','yscale','log')
               
                set(gcf,'PaperPositionMode','manual');
                set(gcf,'PaperUnits','inches');
                set(gcf,'PaperOrientation','portrait');
                set(gcf,'PaperPosition',[.05 .05 7 5]);
                
                end % ifplot
                
            case 3  % =======  COMPLIANCE AND TILT METHOD 3  ======= 
                % rename things
                l1z = lpz;
                l12 = lp2;
                l1p = lp1;
                l2p_1 = l21_p;
                l2z_1 = l2z_p;
                lpz_12 = l1z_p2;
                
                % interpolate transfer functions to the frequencies for the data
                l1z_use = interp1(transfer_fun.f,l1z,fdat);
                l12_use = interp1(transfer_fun.f,l12,fdat);
                l2z_1_use = interp1(transfer_fun.f,l2z_1,fdat);
                l1p_use = interp1(transfer_fun.f,l1p,fdat);
                l2p_1_use = interp1(transfer_fun.f,l2p_1,fdat);
                lpz_12_use = interp1(transfer_fun.f,lpz_12,fdat);

                specZ_1 = spec_Z-((l1z_use)'.*spec_H1); % remove 1 from Z ==> Zcor1
                spec2_1 = spec_H2-(l12_use'.*spec_H1); % remove 1 from 2 ==> 2cor1
                specZ_12 = specZ_1-(l2z_1_use'.*spec2_1); % remove 2cor1 from Zcor1 ==> Zcor12
                specP_1 = spec_P-(l1p_use'.*spec_H1); % remove 1 from P ==> Pcor1
                specP_12 = specP_1-(l2p_1_use'.*spec2_1); % remove 2cor1 from Pcor1 ==> Pcor12
                specZ_12P = specZ_12-(lpz_12_use'.*specP_12); % remove Pcor12 from Zcor12 => Zcor12P

                % power spectra
                poweror = abs(spec_Z).^2*2/(NFFT*dt);
                powerco = abs(specZ_12).^2*2/(NFFT*dt);
                powerc2 = abs(specZ_12P).^2*2/(NFFT*dt);
                
                % transform back to time domain
                ampZ_1 = real(ifft(2*specZ_1,NFFT)./dt);
                ampZ_12 = real(ifft(2*specZ_12,NFFT)./dt);
                ampZ_12P = real(ifft(2*specZ_12P,NFFT)./dt);
                amp_Z_1_nopad = ampZ_1(npad0+1:npad0+npts); %Raw Z data -1
                amp_Z_12_nopad = ampZ_12(npad0+1:npad0+npts); %Raw Z data -12
                amp_Z_12P_nopad = ampZ_12P(npad0+1:npad0+npts); %Raw Z data -12P
                Z_filt_1 = filtfilt(b,a,amp_Z_1_nopad);
                Z_filt_12 = filtfilt(b,a,amp_Z_12_nopad);
                Z_filt_12P = filtfilt(b,a,amp_Z_12P_nopad);
                
                % final corrected data
                Zdat_cor = Z_filt_12P;
                
                %% Some plots
                if ifplot
                % -------- plot progressive corrections to the Z comp ------
                figure(102), clf, set(gcf,'pos',[20 20 1100 700])
                subplot(411)
                title(sprintf('%s  %s',sta,evstr)), hold on
                plot(tt,Z_filt,'-k'); ylabel('Z')
                subplot(412)
                plot(tt,Z_filt_1,'-k'); ylabel('Z-1')
                subplot(413)
                plot(tt,Z_filt_12,'-k'); ylabel('Z-12')
                subplot(414)
                plot(tt,Z_filt_12P,'-k'); ylabel('Z-12P')
                
                set(gcf,'PaperPositionMode','manual');
                set(gcf,'PaperUnits','inches');
                set(gcf,'PaperOrientation','portrait');
                set(gcf,'PaperPosition',[.05 .05 10.5 8]);

                % -------- plot improvements in narrow filter bands ------
                figure(103), clf, set(gcf,'pos',[1050 150 1000 1000])
                for j = 1:length(periods)
                    T1p = periods(j)-0.1*periods(j); T2p= periods(j)+0.1*periods(j);
                    [bp,ap]=butter(2,[1/fNq/T2p,1/fNq/T1p]);
                    Z_filtp  = filtfilt(bp,ap,Zdat);
                    Z_12_filtp = filtfilt(bp,ap,amp_Z_12_nopad);
                    Z_12P_filtp = filtfilt(bp,ap,amp_Z_12P_nopad);

                    subplot(length(periods),3,3*j-2)
                    plot(tt,(Z_filtp),'-k');
                    title(sprintf('orid %s   Z   %d s',strtok(evstr,'-'), periods(j)));

                    subplot(length(periods),3,(3*j)-1)
                    plot(tt,(Z_12_filtp),'-k');
                    title(sprintf('orid %s   Z-12   %d s',strtok(evstr,'-'), periods(j)));

                    subplot(length(periods),3,(3*j))
                    plot(tt,(Z_12P_filtp),'-k');
                    title(sprintf('orid %s   Z-12P   %d s',strtok(evstr,'-'), periods(j)));
                end
                
                set(gcf,'PaperPositionMode','manual');
                set(gcf,'PaperUnits','inches');
                set(gcf,'PaperOrientation','portrait');
                set(gcf,'PaperPosition',[.05 .05 11 12]);
                
                % -------- plot power spectra through correction ------                
                figure(104), clf, hold on
                loglog(fdat,smooth(poweror,40),'-k','LineWidth',.5)
                loglog(fdat,smooth(powerco,40),'-r','LineWidth',.5)
                loglog(fdat,smooth(powerc2,40),'-b','LineWidth',.5)
                legend('Original','RM-tilt','RM-comp')

                title(sprintf('Z-component, Event: %s',evstr));
                xlabel('Frequency (Hz)')
                ylabel('Power Spectrum')
                set(gca,'xscale','log','yscale','log')
                
                set(gcf,'PaperPositionMode','manual');
                set(gcf,'PaperUnits','inches');
                set(gcf,'PaperOrientation','portrait');
                set(gcf,'PaperPosition',[.05 .05 7 5]);
                
                end % ifplot

        end % switch method_trans_func
        

        %% Save figures
        if ifplot && ifsavefigs
            if cor_opt==0, str = 'all'; elseif cor_opt==1, str = 'data'; elseif cor_opt==2, str = 'noise'; end
            save2png_basic(102,[sta,'_',evstr,'_',str,'_Zcor_stages_','meth',num2str(method_trans_func)],[specdir,sta])
            save2png_basic(103,[sta,'_',evstr,'_',str,'_Zcor_filt_','meth',num2str(method_trans_func)],[specdir,sta])
            save2png_basic(104,[sta,'_',evstr,'_',str,'_Zcor_power_','meth',num2str(method_trans_func)],[specdir,sta])
        end
        
        %% Manual QC
        if ifmanualQC
        % >>>>> comment out to not do manual checks
        yn = input('keep this corrected data? [ENTER for yes, otherwise no] ','s');
        if ~isempty(yn) % SKIPPING - but have to make sure no prev. version is carried through
            if isfield(data,'noise_cor') % only bother if a prev. noise_cor has been input
                if ~isempty(data.noise_cor) % only bother if a prev. noise_cor has been saved
                    data.noise_cor = [];
                    data.rmcomp = false;
                    data.rmtilt = false;
                    data = orderfields(data,{'station','network','chans',...
                                         'samprate','nsamps','gcarc','seaz','esaz',...
                                         'phases','dat','tt','raw','noise_cor',...
                                         'NEZ','rmresp','rmtilt','rmcomp'});
                    save(datfile,'data');
                end % if data.noise_cor not empty
            end % if data.noise_cor is a field
            
            if datinfo(is).rmcomp==true || datinfo(is).rmtilt == true
                datinfo(is).rmcomp = false;
                datinfo(is).rmtilt = false;
                save(datinfofile,'datinfo');
            end % if datinfo did show rmcomp or rmtilt;
                
            fprintf(' chose to skip.\n')
            continue
        end
        end
            
        % <<<<< comment out to here
        
        %% Put corrected Z into data structure
        % noise_removal parms
        parms = struct('method',method_trans_func,'specfile',specfile,...
                       'smooth_trans_func',smooth_trans_func,'npt_mav',npt_mav_stf,...
                       'filt_trans_func',filt_trans_func,'low_pass',lp);
        % noise-corrected structure
%         noise_cor = struct('tt',tt,'Zdat',Zdat_cor,'samprate',resamprate,...
%                            'parms',parms,'TF',transfer_fun); % save the TF too! - but data heavy
        noise_cor = struct('tt',tt+evtimes(ie),'Zdat',Zdat_cor,'samprate',resamprate,...
                           'parms',parms);
        
        % put into data structure
        data.noise_cor = noise_cor;
        if method_trans_func == 1, data.rmcomp = true; data.rmtilt = true; end
        if method_trans_func == 2, data.rmcomp = true; end
        if method_trans_func == 3, data.rmcomp = true; data.rmtilt = true; end
        
        % silly thing to get fields in nicer order
        data = orderfields(data,{'station','network','chans',...
                                 'samprate','nsamps','gcarc','seaz','esaz',...
                                 'phases','dat','tt','raw','noise_cor',...
                                 'NEZ','rmresp','rmtilt','rmcomp'});

        %% SAVE
        % datafile
        save(datfile,'data');
        % datinfofile
        if method_trans_func == 1, datinfo(is).rmcomp = true; datinfo(is).rmtilt = true; end
        if method_trans_func == 2, datinfo(is).rmcomp = true; end
        if method_trans_func == 3, datinfo(is).rmcomp = true; datinfo(is).rmtilt = true; end
        save(datinfofile,'datinfo');
        fprintf(' corrected data saved.\n')
   

        
    end % loop on stas

	%% sum up
    fprintf(' STA  CHAN  NEZ  resp  spectra  tilt  comp\n')
    for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f      %1.0f       %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).spectra,datinfo(is).rmtilt,datinfo(is).rmcomp); end
            
end
