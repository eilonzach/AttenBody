% Cycle through events and calculate spectral ratios for all stations.
% 
% This script is a necessary pre-requisite to calculating attenuation using
% the comb of filters (our preferred method) but along the way will
% calculate the attenuation using amplitude spectral ratios (where the
% spectra are calculated using the Thomson multitaper method using MATLAB
% function pmtm). Other useful things to come out of this are the
% signal-to-noise ratio and the frequency at which the signal crosses the
% pre-arrival noise (above which we should not use the signal).
% 
% Z. Eilon 2016
% Updated 02/2019, 07/2019
close all

%% parameters
phase = 'SKS'; % P, S, SKS
component = 'R'; %'Z', 'T', or 'R'
resamprate = 10 ; % new, common sample rate
filtfs = 1./[40 .5]; % [flo fhi] = 1./[Tmax Tmin] in sec
taperx = 0.2;
datwind = [-160 165]; % window of data in eqar structure
sspecwind = [-5 15]; % for P
% sspecwind = [-5 20]; % for S
nspecwind = [-150 -10];
snrmin = 3;
mavwind = 5; % length of moving average to smooth spectrum (1 for no smoothing)
lofrq = 0.05; % uppermost freq to fit (Hz)
% hifrq = 1; % uppermost freq to fit (Hz)
firstev = 815; % SKS got to 1044; S got to 1287; % need to do 160-335 again on all

overwrite = false;
ifplot    = false;
ifsave    = true;

% % project details
% dbname = 'EARdb';
% dbdir = '/Users/zeilon/Work/EastAfrica/EARdb/'; % include final slash

% project details
dbname = 'FRES_PILOT';
dbdir = '~/Dropbox/Work/FRES_PILOT/'; % include final slash



% SNR for spectra
if strcmp(phase(1),'P')
    spec_snr_hi = snrmin; % min SNR for specss/specns sweeping up
    spec_snr_lo = snrmin/2; % min SNR for specss/specns sweeping down
    fsweephi = 0.5; % Hz to start sweeping up looking for crossing
    fsweeplo = 0.2; % Hz to start sweeping down looking for crossing
elseif strcmp(phase(1),'S')
    spec_snr_hi = 1.5; % min SNR for specss/specns sweeping up
    spec_snr_lo = 1.5; % min SNR for specss/specns sweeping down
    fsweephi = 0.17; % Hz to start sweeping up looking for crossing
    fsweeplo = 0.17; % Hz to start sweeping down looking for crossing
end
    
    
%% Preliminaries
wd = pwd;
addpath('matguts')
cd(dbdir);
run([dbdir,dbname,'_startup.m']);
spd = 24*60*60; % seconds per day - to convert time in seconds to serial time

%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

% event details
load([infodir,'/events'],'evinfo'); 
% station details
load([infodir,'/stations'],'stainfo');


for ie = firstev:evinfo.norids %norids % loop on orids
%     if  mags(ie)<6.9, continue, end
    tic
    orid = evinfo.orids(ie);
    fprintf('\n Orid %.0f %s \n\n',orid,evinfo.evtimes_IRISstr{ie}(1:end-4))
    % name files and directories
    evdir = [num2str(orid,'%03d'),'_',evinfo.datestamp(ie,:),'/'];
    datinfofile = [data_eqar_dir,evdir,'_datinfo_',phase];
    arfile      = [data_eqar_dir,evdir,'_EQAR_',phase,'_',component];
      
    % check files exist
    if ~exist([datinfofile,'.mat'],'file'), fprintf('No station mat files for this event\n');continue, end
    if ~exist([arfile,'.mat'],'file'), fprintf('No arfile for this event + phase\n');continue, end
    
    % load files
    load(datinfofile) % loads datinfo stucture
    load(arfile)      % loads eqar structure
    
    % options to skip
    if isempty(datinfo), fprintf('No station mat files for this event\n'); continue, end
    if ~isfield(datinfo,'xcor') 
        fprintf('Need to xcor arrival time for this phase and event\n'), continue
    end
	if ~any([datinfo.xcor]==true)
        fprintf('Need to xcor arrival time for this phase and event\n'), continue
	end
    if isfield(datinfo,'specR')
        if any([datinfo.specR]==true)
            if ~overwrite
                yn = input('delta-tstar already done - overwrite? [y/n] ','s'); 
                if ~strcmp(yn,'y'), fprintf('skipping...\n'), continue, end
            end
        end
    end
    
    %% ------------------ GRAB DATA IN RIGHT FORMAT ------------------
    % prep data structures
    nstas = length(datinfo);
    all_dat0  = zeros(resamprate*diff(datwind),nstas);

    % LOOP ON STAS
    for is = 1:nstas
        fprintf('Station %s... ',datinfo(is).sta)        % APPLY DATA QUALITY CONDITIONS
        if isempty(eqar(is).(['dat' component])),fprintf('no data for this component\n'),continue; end   
        if datinfo(is).xcor==false,fprintf('no xcor for this station+component\n'),continue; end
        if isempty(eqar(is).dT),fprintf('no xcor for this component\n'),continue; end
        if isnan(eqar(is).dT),fprintf('no xcor for this component\n'),continue; end

        % SHIFT TRACES USING XCOR ARRIVAL TIME 
        % shift so arrival is at time=0
        att = (eqar(is).tt-eqar(is).abs_arrT)*spd; % shift to seconds since absolute arrival
        att = round(att,3);
        ja = (att >= datwind(1)) & (att < datwind(2)); % excerpt times according to datwind
        
        % GRAB DESIRED COMPONENT
        switch component
            case 'Z', all_dat0(:,is) = eqar(is).datZ(ja);
            case 'R', all_dat0(:,is) = eqar(is).datR(ja);
            case 'T', all_dat0(:,is) = eqar(is).datT(ja);
        end
        fprintf('got data\n')
    end % loop on stas
    
    % CLEAN DATA for signal and noise
    cp = struct('samprate',resamprate,'pretime',-datwind(1),'prex',-sspecwind(1),'postx',sspecwind(2),...
            'taperx',taperx,'fhi',filtfs(2),'flo',filtfs(1),'npoles',2,'norm',0);
    cp_n = struct('samprate',resamprate,'pretime',-datwind(1),'prex',-nspecwind(1),'postx',nspecwind(2),...
            'taperx',taperx,'fhi',filtfs(2),'flo',filtfs(1),'npoles',2,'norm',0);

    
    [ all_datwf,all_datf,all_datwc,all_datc,~,ttws,tts ] = data_clean( all_dat0,cp );
    [ all_datwf_n,all_datf_n,all_datwc_n,all_datc_n,~,ttws_n,tts_n ] = data_clean( all_dat0,cp_n );
    
    % CALCULATE SNR 
    % from ratio of variance of noise to signal
    snrwf = var(all_datwf)./var(all_datwf_n);
    eqar(1).snr_wf = [];
    eqar = dealto(eqar,'snr_wf',snrwf);

    % ONLY USE GOOD TRACES
    indgd = 1:size(eqar);
    indgd(mean(abs(all_datwf(:,indgd)))==0)     = []; % kill zero traces
    indgd(isnan(mean(abs(all_datwf(:,indgd))))) = []; % kill nan traces
    indgd(snrwf(indgd)<snrmin)                  = []; % kill low snr traces
    indgd(strcmp({eqar(indgd).sta},'DESE'))     = []; % kill DESE - weird station
    if length(indgd) < 2, fprintf('NO GOOD TRACES/ARRIVALS, skip...\n'), continue, end
    
    
    %% ------------------ CALCULATE SPECTRA ------------------
    % WORK OUT SOME KEY VALS
 	wlen_s = diff(sspecwind)*resamprate*(1+4*taperx); % signal window length, accounting for taper+padding
 	wlen_n = diff(nspecwind)*resamprate*(1+4*taperx); % noise window length, accounting for taper+padding
    dt = 1./resamprate;
    ntap=2;
    nft_s=2^nextpow2(wlen_s);
    nft_n=2^nextpow2(wlen_n);
    nyq=0.5*resamprate;

%     hw = waitbar(0,'PROGRESS THROUGH CALCULATION OF SPECTRA');
    for ig = 1:length(indgd)
        is = indgd(ig);
        
        [specn,frqn]=pmtm(all_datwf_n(:,is),ntap,nft_n,resamprate); %specnnn = specn;
        [specs,frq] =pmtm(all_datwf(:,is),  ntap,nft_s,resamprate);
        
        frq=frq(2:length(frq)); % ignore 0 freq
        
        % subset noise spectrum frequencies to those in the signal spectrum
        [~,insf,~] = intersect(frqn,frq);
        specn = moving_average(specn,nft_n/nft_s); % anti-alias for interp spectrum
        specn=specn(insf);
        
        specs=specs(2:end); % ignore 0 freq
         

        eqar(is).frq = frq;
% %         %convert power to amplitude, and integrate to displacement
% %         eqar(is).specn=(specn.^0.5)./(2.*pi.*frq);
% %         eqar(is).specs=(specs.^0.5)./(2.*pi.*frq);
        %convert power to amplitude ALREADY IN DISPLACEMENT
        eqar(is).specn=(specn.^0.5);
        eqar(is).specs=(specs.^0.5);
        % smooth spectra
        eqar(is).specns = moving_average(eqar(is).specn,mavwind);
        eqar(is).specss = moving_average(eqar(is).specs,mavwind); 
        
%         % while testing
%         figure(5); clf; hold on; 
%         plot(frqn,specnnn.^0.5,'--');
%         plot(frq,[eqar(is).specns,eqar(is).specss])
%         set(gca,'yscale','log')
       

        % look for first crossing above fminhi, using SNR = specsnrhi
        eqar(is).fcrosshi = frq(max([find(eqar(is).specss < spec_snr_hi*eqar(is).specns & frq > fsweephi,1,'first'),2])-1); 
        % look for first crossing below fminhi, using SNR = specsnrlo
        eqar(is).fcrosslo = frq(find(eqar(is).specss < spec_snr_lo*eqar(is).specns & frq < fsweeplo,1,'last')); 
        
        
        if isempty(eqar(is).fcrosshi), eqar(is).fcrosshi = 0; end
        if isempty(eqar(is).fcrosslo), eqar(is).fcrosslo = lofrq; end
        
%         % while testing
%         figure(55); clf; set(gcf,'pos',[48 770 1266 497]);hold on
%         semilogy(frq,eqar(is).specn,':r')
%         semilogy(frq,eqar(is).specns,'r-')
%         semilogy(frq,eqar(is).specs,'k:')
%         semilogy(frq,eqar(is).specss,'k-')
%         semilogy(frq,snrmin*eqar(is).specns,'b--')
%         semilogy(eqar(is).fcrosshi*[1 1],max(eqar(is).specss)*[0.001,1],'g');
%         semilogy(eqar(is).fcrosslo*[1 1],max(eqar(is).specss)*[0.001,1],'g');
%         set(gca,'YScale','log','xlim',[0 1.7])
%         title(sprintf('f_{lo} = %.2f     f_{hi} = %.2f',eqar(is).fcrosslo,eqar(is).fcrosshi))
%         
%       pause
%         waitbar(ig/length(indgd),hw)
    end
%     delete(hw)
    
    %% ONLY USE GOOD TRACES
    % find traces with fewer than 4 points in freq space
    bad_f_range = ([eqar(indgd).fcrosshi] - [eqar(indgd).fcrosslo]) < 4*frq(1);
    indgd(bad_f_range)  	= []; % kill traces that have no/little signal above threshold
    if length(indgd) < 3, fprintf('NOT ENOUGH GOOD TRACES/ARRIVALS, skip...\n'), 
        % reset in case populated
        eqar =  dealto(eqar,'dtstar_specR',nan);eqar = dealto(eqar,'std_dtstar_specR',nan);
        [datinfo.dtstar_specR] = deal(false);
        save(arfile,'eqar')
        save(datinfofile,'datinfo')
        continue
    end
    
    %% reset fmax as per crossing frequencies
    hifrq = nanmean([eqar(indgd).fcrosshi],2);
    lofrq = nanmean([eqar(indgd).fcrosslo],2);
    
 	%% ------------------ CALCULATE DT-STAR FROM REFSTA ------------------
%     refsta = 'WISH';
%     iref = find(strcmp({eqar.sta},refsta));
%     refspecss = eqar(iref).specss;
%     
%     lnR_all = zeros(length(eqar(is).specss),length(indgd));
% 	ind = frq<hifrq;
%     figure(78), clf, set(gcf,'position',[440 0 1000 1400]), hold on
%     for ig = 1:length(indgd)
%         is = indgd(ig);
%         ispecss = eqar(is).specss;
%         lnR = log(ispecss./refspecss);
%         fo = fit(frq(ind),lnR(ind),'poly1');
%         
%         subplot(2,1,1), hold on
%         hr = plot(frq,lnR,'Linewidth',1.5);
%         hrf = plot(fo);
%         xlim([0 hifrq])
%         
%         subplot(2,1,2), hold on
%         hs = plot(ttws,all_datwf(:,is)./max(max(abs(all_datwf(:,indgd))))+eqar(is).gcarc,'Linewidth',1.5);
%         
%         set(hs,'color',colour_get(eqar(is).slon,max([eqar.slon]),min([eqar.slon])))
%         set(hr,'color',colour_get(eqar(is).slon,max([eqar.slon]),min([eqar.slon])))
%         set(hrf,'color',colour_get(eqar(is).slon,max([eqar.slon]),min([eqar.slon])))
%         
%         lnR_all(:,ig) = lnR;
%     end
    
    %% ------------------ CALCULATE DIFFERENTIAL T-STAR ------------------
    specss = zeros(length(eqar(is).specss),length(indgd));
    for ig = 1:length(indgd)
        is = indgd(ig);
        specss(:,ig) = eqar(is).specss;
    end
    
    % CALC DELTA-TSTAR.
    fprintf('Calculate least-squares differential t-star\n')
    [ delta_tstar,cov_dtstar,std_dtstar ] = xspecratio( specss,frq,[eqar(indgd).fcrosshi],[eqar(indgd).fcrosslo],1,ifplot );
    
    
    %% ---------------------- STORE RESULTS -----------------------
    % STORE RESULTS
    fprintf('Recording results in arrival structure...')
    % prep eqar to receive new fields
    eqar(1).dtstar_specR = []; eqar(1).std_dtstar_specR = []; eqar(1).par_dtstar_specR = [];
    par_dtstar = struct('comp',component,'filtfs',filtfs,'swindow',sspecwind,'nwindow',nspecwind,...
                          'taperx',taperx,'mavwind',mavwind,...
                          'spec_snr_hi',spec_snr_hi,'spec_snr_lo',spec_snr_lo,'snrmin',snrmin);
                      
    if isfield(eqar,'dtstar') % relabel from old convention
        for isdt = 1:length(eqar), eqar(isdt).dtstar_specR = eqar(isdt).dtstar; end; eqar = rmfield(eqar,'dtstar');
        for isdt = 1:4, eqar(isdt).std_dtstar_specR = eqar(isdt).std_dtstar; end; eqar = rmfield(eqar,'std_dtstar');
    end

    eqar(indgd) =  dealto(eqar(indgd),'dtstar_specR',delta_tstar);
    eqar(indgd) =  dealto(eqar(indgd),'std_dtstar_specR',std_dtstar);
    eqar(indgd) =  dealto(eqar(indgd),'par_dtstar_specR',par_dtstar);

	indbd = setdiff(1:length(eqar),indgd);
    if ~isempty(indbd)
    eqar(indbd) =  dealto(eqar(indbd),'dtstar_specR',nan);
    end
    
    if any(abs(delta_tstar)>5)
        weird
    end
    
    
    %% -------------------------- PLOTS ---------------------------
    if ifplot
    %% look at maximum frequencies above noise
    figure(87), clf, hold on
    fcross = [eqar(indgd).fcrosshi]';
    hist(fcross,100); xlim([0 2])
    line(hifrq*[1 1],[0 max(hist(fcross))],'Color','r','LineStyle','--','LineWidth',2)
    title('maximum freq where signal is above noise')

    %% spectral plot
    figure(150), clf, set(gcf,'position',[16 132 898 1213]), hold on
    for ig = 1:length(indgd)
        is = indgd(ig);

        % plot data series
        subplot(212), hold on
        nrmamp = max(abs([all_datwf_n(:,is);all_datwf(:,is)]));
        hpd = plot(ttws_n,all_datwf_n(:,is)./nrmamp,ttws,all_datwf(:,is)./nrmamp,'linewidth',2);
%         set(hpd,'color',colour_get(eqar(is).slon,max([eqar.slon]),min([eqar.slon])))
        set(hpd,'color',colour_get(eqar(is).dtstar_specR,1,-1))
        xlim([-40 sspecwind(2)]);
        xlabel('time from pick, s')

        % plot spectrum
        subplot(211), hold on
%         n4=length(find(frq<=0.8*nyq));
        usefrq = find(frq<=eqar(is).fcrosshi & frq>=eqar(is).fcrosslo);
        hps = semilogx(frq(usefrq),log10(eqar(is).specss(usefrq)),'k','LineWidth',1); % << here using smoothed spectrum 
        hpn = semilogx(frq(usefrq),log10(eqar(is).specns(usefrq)),'k','LineWidth',1.5); % << here using smoothed spectrum 
%         set(hps,'color',colour_get(eqar(is).slon,max([eqar.slon]),min([eqar.slon])))
        set(hpn,'color',colour_get(eqar(is).dtstar_specR,1,-1))
%         set(hps,'color',colour_get(eqar(is).dtstar,1,-1,jet))
        % xlim(ax(1),[0 0.8*nyq]);
        xlabel('Hz'); 
        ylabel('amplitude, nm/Hz')
        % a=axis(ax(1));
        % % ylim_max=10^(floor(log10(max(specs(1:n4))))+1);
        % % ylim(ax(1),ylim_max.*[1.e-4 1])
        xlim([0 hifrq]);
%         ylim([-9,-2]);
% pause
    end % loop on good stas
    
    %% quick and dirty geog plot
%     figure(4)
%     scatter([eqar(indgd).slon]',[eqar(indgd).slat]',200,[eqar(indgd).dtstar]','filled')
%     figure(5)
%     plot([eqar(indgd).fcrosshi]',[eqar(indgd).dtstar]','o')

    eval(sprintf('plot_ATTEN_TandF_domain_%s( eqar(indgd))',dbname ))
    end % ifplot
    
	%% -------------------------- SAVE ---------------------------
    if ifsave
    % SAVE
    save(arfile,'eqar')
    fprintf(' saved\n')
    % RECORD DTSTARCALC IN DATINFO 
    [datinfo.specR] = deal(false);
	[datinfo(indgd).specR] = deal(true);
    save(datinfofile,'datinfo')
    
    fprintf(' STA  CHAN  NEZ  resp  tilt  comp  xcor  specR\n')
	for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f     %1.0f     %1.0f     %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).rmtilt,datinfo(is).rmcomp,datinfo(is).xcor,datinfo(is).specR); end
    end

    
    toc  
end % loop on orids

% results_PARSE