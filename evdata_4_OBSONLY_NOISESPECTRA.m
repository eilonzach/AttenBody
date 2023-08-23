% Script to compute spectral characteristics of noise before each event.
% This script goes through data event by event, then station by station.
% For OBS stations, it pulls out the 12 hours of noise prior to each event
% and computes spectra on a set of Nwin windows of noise, saving the
% outputs into structures within the data directories
clear all

overwrite = false;

% day info to calculate tilt and compliance
noisewind = 43200;     % length of noise window before event time (sec).
T    = 6000;  % the length of each time window, in sec  
Nwin = 10;    % N of time windows into which to chop the noise 

resamprate = 5; % NEW SAMPLE RATE TO DOWNSAMP TO

ifQC = 1;
ifplot = 1; 
ifsavefigs = 1;


%% paths
cd('/Users/zeilon/Documents/MATLAB/CASC_atten/')
addpath('/Users/zeilon/Documents/MATLAB/helen_dataprocessing')
addpath('matguts')
% path to sensor tilt data
tiltdir = '~/Work/CASCADIA/CAdb/tilt_compliance/';
% antelope db details
dbdir = '/Users/zeilon/Work/CASCADIA/CAdb/'; % needs final slash
dbnam = 'cascBIGdb';
% path to top level of directory tree for data
datadir = '/Volumes/DATA_mini2/CASCADIA/DATA/'; % needs final slash


%% get to work
if ~(T > noisewind/Nwin), error('Make T larger so there is overlap between windows'); end
[ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam ); % load events data
[ nstas,stas,slats,slons,selevs,ondate,offdate,staname,statype ] = db_stadata( dbdir,dbnam ); %load stations data

for ie = 1:269 % 1:norids % loop on orids
    fprintf('\n Orid %.0f %s \n\n',orids(ie),epoch2str(evtimes(ie),'%Y-%m-%d %H:%M:%S'))
    evdir = [num2str(orids(ie),'%03d'),'_',epoch2str(evtimes(ie),'%Y%m%d%H%M'),'/'];
    evday = epoch2str(evtimes(ie),'%Y%j');
    datinfofile = [datadir,evdir,'_datinfo'];
           
    if any((evtimes(ie)-evtimes)>0 & (evtimes(ie)-evtimes) < 20*60*60)
        fprintf('Another event within prev 20 hrs... skipping\n')
        continue
    end

    % check files exist
    if ~exist([datinfofile,'.mat'],'file'), fprintf('No data at all for this event\n');continue, end
    load(datinfofile)
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end

    for is = 1:length(datinfo) % loop on stas
        sta = datinfo(is).sta; % sta name
        if strcmp(statype(strcmp(stas,sta)),'OBS') ~= 1, continue, end % skip if not an OBS
        
        ofile = [datadir,evdir,sta,'_spectra.mat'];
        datfile = [datadir,evdir,sta,'.mat'];
        
        fprintf('Station %.0f, %s...',is,sta)
        
        %% load data, be v. careful about overwriting etc., w/ redundancy built in
        % conditions before even loading - only do if response removed
        if datinfo(is).spectra & ~overwrite  % skip if already calculated spectra
            fprintf(' done already...')
%             if ifplot
%                 fprintf(' plotting')
%                 plot_1evt1sta_spectra(ofile,ifsavefigs,[dbdir,'spectra/',sta])
%             end
            fprintf('\n'), continue% skip if already calculated spectra
        end
        if datinfo(is).rmresp ~= 1, fprintf(' must remove response\n'), continue, end % skip if response not already removed
        
        % load data
        load(datfile); % load sta data for this evt
        % second pass on conditions
        if data.rmtilt, fprintf(' done already\n'), continue, end % skip if already removed tilt
        if ~data.rmresp, fprintf(' must remove response\n'), continue, end % skip if response not already removed
        
        % outfile + check if already exists
        if exist(ofile,'file')==2
            fprintf(' specfile exists... ')
            if overwrite == false, fprintf('skipping...\n'),continue
            else fprintf('overwriting... '), delete(ofile); end
        end
        
        %% sort out channels        
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

        %% Get data, reseamp
        inds = noisewind*data.samprate;
        
        noiseZ  = data.dat(1:inds,icz); noiseZ(isnan(noiseZ))=0;
        noiseH1 = data.dat(1:inds,icn); noiseH1(isnan(noiseH1))=0;
        noiseH2 = data.dat(1:inds,ice); noiseH2(isnan(noiseH2))=0;
        noiseP  = data.dat(1:inds,ich); noiseP(isnan(noiseP))=0;
        noiset = data.tt(1:inds,1);
        
        % NEED TO USE THE RAW PRESSURE IF THERE IS A LAMONT APG
        if strcmp(which_OBS(sta),'LDEO')
            fprintf('LDEO APG - using raw data ')
        	ichr = find(strcmp(raw.chans.component,'H'));
            noiseP  = data.raw.dat(1:inds,ichr); noiseP(isnan(noiseP))=0;
        end
                    
        % resamp time vector
        tt0 = 0;
        tt1 = noisewind - 1./resamprate;
        tt = [tt0:(1./resamprate):tt1]';
                
        noiseZ  = interp1(noiset-noiset(1),noiseZ,tt); noiseZ(isnan(noiseZ)) = 0;
        noiseH1 = interp1(noiset-noiset(1),noiseH1,tt); noiseH1(isnan(noiseH1)) = 0;
        noiseH2 = interp1(noiset-noiset(1),noiseH2,tt); noiseH2(isnan(noiseH2)) = 0;
        noiseP  = interp1(noiset-noiset(1),noiseP,tt); noiseP(isnan(noiseP)) = 0;
        
        %% plot the spectra of the whole thing
        if ifplot
        [ whole_spectra_fig ] = plot_quick_spectra( noiseZ,noiseH1,noiseH2,noiseP,tt(1:T*resamprate),Nwin);
        pause(.1)
        end 

        %% QC the noise windows
        if ifQC
            fprintf('QC wins... ')
            [ goodwins ] = QC_staspectra_windows( noiseZ,noiseH1,noiseH2,noiseP,tt(1:T*resamprate),Nwin,1);
            fprintf('killed %.0f... ',Nwin-length(goodwins))
            if ifplot && (Nwin-length(goodwins) > 0)
                plot_quick_spectra( noiseZ,noiseH1,noiseH2,noiseP,tt(1:T*resamprate),Nwin,goodwins);
                pause(.1)
            end
        else
            goodwins = 1:Nwin;
        end
        
        %% calc some basics
        oldsamprate = data.samprate;
        samprate = resamprate;
        dt = 1./samprate;
        npts = T*samprate;
        Ppower = nextpow2(npts);
        NFFT = 2^Ppower;
        npad0 = (NFFT-npts)/2;
        f = samprate/2*linspace(0,1,NFFT/2+1)';
        win_t0 = round(linspace(1,noisewind-T+1,Nwin));
        win_pt_start = (win_t0-1)*samprate+1;
        win_pt_end = win_pt_start+npts-1;
        taperf = (win_pt_end(1)-win_pt_start(2))/npts;

        %% scales for seismometer and pressure chans...?
        scalez = 1;
        scalep = 1;

        %% Prep output structures
        % ----- open zeros array for spectrum and conj(spectrum) for 4 channels -----
        spectrum_Z=zeros(Nwin,length(f));
        spectrum_H1=zeros(Nwin,length(f));
        spectrum_H2=zeros(Nwin,length(f));
        spectrum_P=zeros(Nwin,length(f));
        cspectrum_Z=zeros(Nwin,length(f));
        cspectrum_H1=zeros(Nwin,length(f));
        cspectrum_H2=zeros(Nwin,length(f));
        cspectrum_P=zeros(Nwin,length(f));

        c11_stack=zeros(length(f),1);
        c22_stack=zeros(length(f),1);
        cpp_stack=zeros(length(f),1);
        czz_stack=zeros(length(f),1);

        c1z_stack=zeros(length(f),1);
        c2z_stack=zeros(length(f),1);
        cpz_stack=zeros(length(f),1);
        c12_stack=zeros(length(f),1);
        c1p_stack=zeros(length(f),1);
        c2p_stack=zeros(length(f),1);

        C1z_stack = zeros(length(f),1);
        C2z_stack = zeros(length(f),1);
        Cpz_stack = zeros(length(f),1);
        C12_stack = zeros(length(f),1);
        C2p_stack = zeros(length(f),1);
        C1p_stack = zeros(length(f),1);

        Q1z_stack = zeros(length(f),1);
        Q2z_stack = zeros(length(f),1);
        Qpz_stack = zeros(length(f),1);
        Q12_stack = zeros(length(f),1);
        Q1p_stack = zeros(length(f),1);
        Q2p_stack = zeros(length(f),1);

        Nwin_stack = 0;
        %% loop over windows
        for iwin = 1:Nwin % Window index
            if ~ismember(goodwins,iwin) % 
                continue
            end
            j1 = win_pt_start(iwin);
            j2 = win_pt_end(iwin);

       
            amp_Z  = noiseZ(j1:j2);
            amp_Z  = amp_Z.*flat_hanning(tt(j1:j2),taperf*dt*npts);
            amp_Z  = detrend(amp_Z,0);
            amp_Z  = padarray(amp_Z,[npad0 0],'both');
            spectrum = fft(amp_Z,NFFT).*dt;
            spec_Z = spectrum(1:NFFT/2+1);
            cspec_Z = conj(spec_Z);
            spectrum_Z(iwin,1:length(f)) = spec_Z;
            cspectrum_Z(iwin,1:length(f)) = cspec_Z;

            amp_H1  = noiseH1(j1:j2);
            amp_H1  = amp_H1.*flat_hanning(tt(j1:j2),taperf*dt*npts);
            amp_H1  = detrend(amp_H1,0);
            amp_H1  = padarray(amp_H1,[npad0 0],'both');
            spectrum = fft(amp_H1,NFFT).*dt;
            spec_H1 = spectrum(1:NFFT/2+1);
            cspec_H1 = conj(spec_H1);
            spectrum_H1(iwin,1:length(f)) = spec_H1;
            cspectrum_H1(iwin,1:length(f)) = cspec_H1;

            amp_H2  = noiseH2(j1:j2);
            amp_H2  = amp_H2.*flat_hanning(tt(j1:j2),taperf*dt*npts);
            amp_H2  = detrend(amp_H2,0);
            amp_H2  = padarray(amp_H2,[npad0 0],'both');
            spectrum = fft(amp_H2,NFFT).*dt;
            spec_H2 = spectrum(1:NFFT/2+1);
            cspec_H2 = conj(spec_H2);
            spectrum_H2(iwin,1:length(f)) = spec_H2;
            cspectrum_H2(iwin,1:length(f)) = cspec_H2;

            amp_P  = noiseP(j1:j2);
            amp_P  = amp_P.*flat_hanning(tt(j1:j2),taperf*dt*npts);
            amp_P  = detrend(amp_P,0);
            amp_P  = padarray(amp_P,[npad0 0],'both');
            spectrum = fft(amp_P,NFFT).*dt;
            spec_P = spectrum(1:NFFT/2+1);
            cspec_P = conj(spec_P);
            spectrum_P(iwin,1:length(f)) = spec_P;
            cspectrum_P(iwin,1:length(f)) = cspec_P;

            % power spectrum for each segment
            c11 = abs(spec_H1).^2*2/(NFFT*dt);
            c22 = abs(spec_H2).^2*2/(NFFT*dt);
            cpp = abs(spec_P).^2*2/(NFFT*dt);
            czz = abs(spec_Z).^2*2/(NFFT*dt);
            % complex cross-spectrum for each segment
            c1z = spec_H1.*cspec_Z*2/(NFFT*dt);
            c2z = spec_H2.*cspec_Z*2/(NFFT*dt);
            cpz = spec_P.*cspec_Z*2/(NFFT*dt);
            c12 = spec_H1.*cspec_H2*2/(NFFT*dt);
            c1p = spec_H1.*cspec_P*2/(NFFT*dt);
            c2p = spec_H2.*cspec_P*2/(NFFT*dt);

            % other cross spectral density functions
            C1z = (real(spec_H1).*real(spec_Z)+imag(spec_H1).*imag(spec_Z));
            C2z = (real(spec_H2).*real(spec_Z)+imag(spec_H2).*imag(spec_Z));
            Cpz = (real(spec_P).*real(spec_Z)+imag(spec_P).*imag(spec_Z));
            C12 = (real(spec_H1).*real(spec_H2)+imag(spec_H1).*imag(spec_H2));
            C1p = (real(spec_H1).*real(spec_P)+imag(spec_H1).*imag(spec_P));
            C2p = (real(spec_H2).*real(spec_P)+imag(spec_H2).*imag(spec_P));

            Q1z = (real(spec_H1).*imag(spec_Z)+imag(spec_H1).*real(spec_Z));
            Q2z = (real(spec_H2).*imag(spec_Z)+imag(spec_H2).*real(spec_Z));
            Qpz = (real(spec_P).*imag(spec_Z)+imag(spec_P).*real(spec_Z));
            Q12 = (real(spec_H1).*imag(spec_H2)+imag(spec_H1).*real(spec_H2));
            Q1p = (real(spec_H1).*imag(spec_P)+imag(spec_H1).*real(spec_P));
            Q2p = (real(spec_H2).*imag(spec_P)+imag(spec_H2).*real(spec_P));

            cc=interp1(1:64,colormap('jet'),(iwin/(Nwin+1))*63+1);

            % stack
            %         ind_amp=find((f>fmin).*(f<fmax));
            %         if(mean(c11(ind_amp))<chorz_max)&&((mean(c22(ind_amp))<chorz_max))
            c11_stack=c11_stack+c11;
            c22_stack=c22_stack+c22;
            cpp_stack=cpp_stack+cpp;
            czz_stack=czz_stack+czz;

            c1z_stack=c1z_stack+c1z;
            c2z_stack=c2z_stack+c2z;
            cpz_stack=cpz_stack+cpz;
            c12_stack=c12_stack+c12;
            c1p_stack=c1p_stack+c1p;
            c2p_stack=c2p_stack+c2p;

            C1z_stack = C1z_stack+C1z;
            C2z_stack = C2z_stack+C2z;
            Cpz_stack = Cpz_stack+Cpz;
            C12_stack = C12_stack+C12;
            C1p_stack = C1p_stack+C1p;
            C2p_stack = C2p_stack+C2p;

            Q1z_stack = Q1z_stack+Q1z;
            Q2z_stack = Q2z_stack+Q2z;
            Qpz_stack = Qpz_stack+Qpz;
            Q12_stack = Q12_stack+Q12;
            Q2p_stack = Q2p_stack+Q2p;
            Q1p_stack = Q1p_stack+Q1p;

            Nwin_stack=Nwin_stack+1;
            vec_win(iwin)=1;
            %         end
        end
        
        % delete bad windows
        spectrum_Z = spectrum_Z(goodwins,:);
        spectrum_H1=spectrum_H1(goodwins,:);
        spectrum_H2=spectrum_H2(goodwins,:);
        spectrum_P = spectrum_P(goodwins,:);
        cspectrum_Z = cspectrum_Z(goodwins,:);
        cspectrum_H1=cspectrum_H1(goodwins,:);
        cspectrum_H2=cspectrum_H2(goodwins,:);
        cspectrum_P = cspectrum_P(goodwins,:);
        win_pt_start = win_pt_start(goodwins);
        win_pt_end = win_pt_end(goodwins);

        %% stack over the windows
        %Normalization

        c11_stack=c11_stack/Nwin_stack;
        c22_stack=c22_stack/Nwin_stack;
        cpp_stack=cpp_stack/Nwin_stack;
        czz_stack=czz_stack/Nwin_stack;

        c1z_stack=c1z_stack/Nwin_stack;
        c2z_stack=c2z_stack/Nwin_stack;
        cpz_stack=cpz_stack/Nwin_stack;
        c12_stack=c12_stack/Nwin_stack;
        c1p_stack=c1p_stack/Nwin_stack;
        c2p_stack=c2p_stack/Nwin_stack;

        C1z_stack=C1z_stack/Nwin_stack;
        C2z_stack=C2z_stack/Nwin_stack;
        Cpz_stack=Cpz_stack/Nwin_stack;
        C12_stack=C12_stack/Nwin_stack;
        C1p_stack=C1p_stack/Nwin_stack;
        C2p_stack=C2p_stack/Nwin_stack;

        Q1z_stack=Q1z_stack/Nwin_stack;
        Q2z_stack=Q2z_stack/Nwin_stack;
        Qpz_stack=Qpz_stack/Nwin_stack;
        Q12_stack=Q12_stack/Nwin_stack;
        Q1p_stack=Q1p_stack/Nwin_stack;
        Q2p_stack=Q2p_stack/Nwin_stack;

        % Coherence

        coh1z_stack = abs(c1z_stack).^2./(c11_stack.*czz_stack);
        coh2z_stack = abs(c2z_stack).^2./(c22_stack.*czz_stack);
        cohpz_stack = abs(cpz_stack).^2./(cpp_stack.*czz_stack);
        coh1p_stack = abs(c1p_stack).^2./(c11_stack.*cpp_stack);
        coh2p_stack = abs(c2p_stack).^2./(c22_stack.*cpp_stack);
        coh12_stack = abs(c12_stack).^2./(c11_stack.*c22_stack);

        % Phase

        ph1z_stack = 180/pi.*atan(Q1z_stack./C1z_stack);
        ph2z_stack = 180/pi.*atan(Q2z_stack./C2z_stack);
        phpz_stack = 180/pi.*atan(Qpz_stack./Cpz_stack);
        ph1p_stack = 180/pi.*atan(Q1p_stack./C1p_stack);
        ph2p_stack = 180/pi.*atan(Q2p_stack./C2p_stack);
        ph12_stack = 180/pi.*atan(Q12_stack./C12_stack);

        % Admittance
        
        ad1z_stack = abs(c1z_stack)./c11_stack;
        ad2z_stack = abs(c2z_stack)./c22_stack;
        adpz_stack = abs(cpz_stack)./cpp_stack;
        ad1p_stack = abs(c1p_stack)./c11_stack;
        ad2p_stack = abs(c2p_stack)./c22_stack;
        ad12_stack = abs(c12_stack)./c11_stack;


        %% SAVE
        save(ofile,'sta','evdir','f','tt','dt','NFFT','win_pt_start','win_pt_end','taperf','goodwins',...
            'c11_stack','c22_stack','czz_stack','cpp_stack',...
            'c12_stack','c1z_stack','c2z_stack','c1p_stack','c2p_stack','cpz_stack',...
            'C12_stack','C1z_stack','C2z_stack','C1p_stack','C2p_stack','Cpz_stack',...
            'Q12_stack','Q1z_stack','Q2z_stack','Q1p_stack','Q2p_stack','Qpz_stack',...
            'coh12_stack','coh1z_stack','coh2z_stack','coh1p_stack','coh2p_stack','cohpz_stack',...
            'ph12_stack','ph1z_stack','ph2z_stack','ph1p_stack','ph2p_stack','phpz_stack',...
            'ad12_stack','ad1z_stack','ad2z_stack','ad1p_stack','ad2p_stack','adpz_stack',...
            'spectrum_P','spectrum_Z','spectrum_H1','spectrum_H2',...
            'cspectrum_P','cspectrum_Z','cspectrum_H1','cspectrum_H2');
        fprintf('spectral details saved.\n')
          
        
        %% plot
        if ifplot
        plot_1evt1sta_spectra(ofile,ifsavefigs,[dbdir,'spectra/',sta])
%         pause
        end
        %% log spectra calc and save

        datinfo(is).spectra = true;
        save(datinfofile,'datinfo')
    end % loop on stas

	%% sum up
    fprintf(' STA  CHAN  NEZ  resp  spectra  tilt  comp\n')
    for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f      %1.0f       %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).spectra,datinfo(is).rmtilt,datinfo(is).rmcomp); end
            
end
