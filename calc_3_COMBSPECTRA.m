% cycle through events and calculate attenuation for all stations by
% running a comb of filters over each pair of stations' traces and
% calculating relative phase and differential amplitude at each frequency,
% before fitting the spectra simultaneously with curves whose slopes relate
% simply to delta-tstar
close all

% project details
dbname = 'EARdb';
dbdir = '/Users/zeilon/Dropbox/Work/EARdb/'; % include final slash

% project details
% dbname = 'FRES_PILOT';
% dbdir = '~/Dropbox/Work/FRES_PILOT/'; % include final slash

% dbname = 'AKTAdb';
% dbdir = '/Users/zeilon/Work/AKTA_Cristhian/AKTAdb/'; % include final slash


%% parameters
phase = 'P';
component = 'Z'; %'Z', 'R', or 'T'

overwrite  = true; % should be "true" if ifplotonly 
ifsave     = true;
ifplot     = false; % 2 for the individual combs...
ifplotonly = false;

% datproc parms
resamprate = 4 ; % new, common sample rate
filtfs = 1./[40 1]; % [flo fhi] = 1./[Tmax Tmin] in sec
taperx = 0.2;
datwind = [-160 165]; % window of data in eqar structure
specwind = [-15 15]; % [-prex postx] So both values are relative to arrival (prex is thus +ive)
snrmin = 4;

firstorid = 2299; % 412 is a good one for the amp2phiwt working, 2200

% combspec parms
parms.comb.Tmin = 1;
parms.comb.Tmax = 25;
parms.comb.Nwds = 20;
parms.comb.Tw_opt = 'scale';
parms.comb.npol = 4;
parms.comb.resamprate = resamprate;

parms.wind.pretime = -datwind(1); % seconds before the predicted arrival that the data time series starts. Should be a positive number
parms.wind.prex = -specwind(1); % seconds before the predicted arrival to take window. So if you want to start before the phase arrive, this should be positive. 
parms.wind.postx = specwind(2); % seconds after the predicted arrival to take window. So if you want to start after the phase arrive, this should be positive. 
parms.wind.taperx = 0.1;

parms.qc.minacor = 0.6;
parms.qc.maxphi = 5;
parms.qc.mtmRcomp = true; % option to compare comb and mtm estimates of As
parms.qc.focus_threshold = 3; % threshold value for focus metric to flag up as bad

parms.inv.amp2phiwt = 2; % ratio of upweight of amplitude vs. phase spectra fit

parms.inv.fmin_cskip = 0.1; 
parms.inv.fmin = 0.045; % caps maximum value otherwise set by fcrosslo
parms.inv.fmax = 1;    % caps minimum value otherwise set by fcrosshi
parms.inv.corr_c_skip = true; % correct for cycle-skip. 
parms.inv.ifwt = true; % weight data by xcorr value in fit of narrow-band filtered data between two stations
parms.inv.R2default = 0.9; % if zero, do not weight by R2 of each pairwise fit, Else, weight by multiple of this default.
parms.inv.opt = 1; %USE 1!   for calc_fdependent...  1 is all in one method; 2 is one-by-one method

parms.inv.alpha = 0.27; % assumed frequency dependency of attenuation - set to empty to test for this

if strcmp(dbname,'AKTAdb')
    datwind = [-180 90]; % window of data in eqar structure
    resamprate = 5 ; % new, common sample rate
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


for ie = firstorid:evinfo.norids %evinfo.norids % 335:norids % loop on orids
%     if  mags(ie)<6.9, continue, end
    tic
    orid = evinfo.orids(ie);
    fprintf('\n Orid %.0f %s %s-%s\n\n',orid,evinfo.evtimes_IRISstr{ie}(1:end-4),phase,component)
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

    % bring "dtstar" vs. "specR" notation into alignment
    if isfield(datinfo,'dtstar') && ~isfield(datinfo,'specR')
        [datinfo.specR] = deal(datinfo.dtstar);
        datinfo = rmfield(datinfo,'dtstar');
    end

    % options to skip
    if isempty(datinfo), fprintf('No station mat files for this event\n'); continue, end
    if ~isfield(datinfo,'xcor') 
        fprintf('Need to xcor arrival time for this phase and event\n'), continue
    elseif ~any([datinfo.xcor]==true)
        fprintf('Need to xcor arrival time for this phase and event\n'), continue
	end
    if ~isfield(datinfo,'specR') && ~isfield(datinfo,'dtstar')
        fprintf('Need to calc specR for this phase and event\n'), continue
    else
        if (isfield(datinfo,'specR') && ~any([datinfo.specR]==true)) || ...
           (isfield(datinfo,'dtstar') && ~any([datinfo.dtstar]==true))
            fprintf('Need to calc specR for this phase and event\n'), continue
        end
    end
    if isfield(datinfo,'comb')
        if any([datinfo.specR]==true)
            if ~ifplotonly
            if ~overwrite 
                yn = input('delta-tstar_comb already done - overwrite? [y/n] ','s'); 
                if ~strcmp(yn,'y'), fprintf('skipping...\n'), continue, end
            end
            end
        end
    end
    
    if ifplotonly
        indgd = ~isnan([eqar.dtstar_comb]);
        if ~any(indgd), continue; end
%         if ~any([eqar(indgd).slat]>12.5 & [eqar(indgd).slon]<42  & [eqar(indgd).slon]>39), continue; end
        plot_ATTEN_TandF_domain_COMB_EAR(eqar)
        pause
        continue
    end

    % reset all comb-related fields in eqar if saving
    if ifsave
        for fnames = {'dtstar_comb','dT_comb','A0_comb','alpha_comb'}
            eqar = dealto(eqar,fnames{:},NaN);
        end
    end

    %% ------------------ GRAB DATA IN RIGHT FORMAT ------------------
    % prep data structures
    nstas = length(datinfo);
    all_dat0  = zeros(unique([eqar.samprate])*diff(datwind),nstas);
    
    % calc dc timeshift from abs to pred
    yx = false(nstas,1);
    for is = 1:nstas, yx(is)=~isempty(eqar(is).abs_arrT); end
    dcTshft = nanmean([eqar(yx).pred_arrT] - [eqar(yx).abs_arrT])*spd;

    % LOOP ON STAS
    for is = 1:nstas
        fprintf('Station %s... ',datinfo(is).sta)        % APPLY DATA QUALITY CONDITIONS
        if isempty(eqar(is).(['dat' component])),fprintf('no data for this component\n'),continue; end   
        if isempty(eqar(is).dT),fprintf('no xcor for this component\n'),continue; end
        if isnan(eqar(is).dT),fprintf('no xcor for this component\n'),continue; end
        try
            if isnan(eqar(is).dtstar_specR),fprintf('no specR-dtstar for this component\n'),continue; end
        catch
            warning('eqar seems to call specR method "dtstar"')
            if isnan(eqar(is).dtstar),fprintf('no specR-dtstar for this component\n'),continue; end
        end
        if strcmp(eqar(is).sta,'M08C'), eqar(is).slon = -124.895400; end

        % SHIFT TRACES USING Predicted ARRIVAL TIME 
        att = (eqar(is).tt-eqar(is).pred_arrT)*spd + dcTshft; % shift to since predicted arrival 
        % shift so arrival is at time=0
        ja = (att >= datwind(1)) & (att < datwind(2)); % excerpt times according to datwind
        
        % GRAB DESIRED COMPONENT
        switch component
            case 'Z', all_dat0(:,is) = eqar(is).datZ(ja);
            case 'R', all_dat0(:,is) = eqar(is).datR(ja);
            case 'T', all_dat0(:,is) = eqar(is).datT(ja);
        end
        fprintf('got data\n')
    end % loop on stas

    % find OBS stas
    isob = zeros(size(eqar)); for is = 1:length(eqar), isob(is) = ~isempty(which_OBS(eqar(is).sta)); end
    
    % ONLY USE GOOD TRACES
    indgd = 1:length(eqar);
    indgd(mean(abs(all_dat0(:,indgd)))==0)     = []; % kill zero traces
    indgd(mean(abs(all_dat0(:,indgd)))<1e-12)     = []; % kill zero traces
    indgd(isnan(mean(abs(all_dat0(:,indgd))))) = []; % kill nan traces
    if ~isfield(eqar,'snr_wf'), if overwrite,keyboard; end, continue; end % skip if snr not even calculated
    indgd([eqar(indgd).snr_wf]<snrmin)                  = []; % kill low snr traces
    indgd([eqar(indgd).acor] < median([eqar(indgd).acor]) - 2*std([eqar(indgd).acor])) = []; % kill traces with acor that (somehow) passed QC but is actually more than 2std lower than median
    if length(indgd) < 2, fprintf('NO GOOD TRACES/ARRIVALS, skip...\n'), 
            if ifsave && overwrite
                fprintf('Wiping all comb and saving eqar\n')
                save(arfile,'eqar')
            end
        continue
    end
    

    % only keep good stas
    all_dat0_gd = all_dat0(:,indgd);

    % get spectral ratios for comparison
    mtmspecss = struct('frq',eqar(indgd(1)).frq,'specs',[eqar(indgd).specs]);
    
    %% GET station-by-station acceptable f range from crossing freqs
    fcross = [[eqar(indgd).fcrosslo]',[eqar(indgd).fcrosshi]'];
    fcross(fcross(:,1) < parms.inv.fmin,1) = parms.inv.fmin;
    fcross(fcross(:,2) > parms.inv.fmax,2) = parms.inv.fmax;
    
    %% resamp - speeds up the combing
    fprintf('resampling to %.0f Hz\n',resamprate);
    if 1/resamprate>0.5*parms.comb.Tmin, error('resamprate too small for the stated highest comb-filter\n'); end
    tt0 = att(ja)';
    tt1 = [tt0(1):1/resamprate:tt0(end)]';
    nsamps = length(tt1); nstas = size(all_dat0_gd,2);
    all_dat0_resamp = zeros(nsamps,nstas);
    for is = 1:nstas
        all_dat0_resamp(:,is) = interp1(tt0,all_dat0_gd(:,is),tt1);
    end

    %% ROUND 1 - RUN THROUGH COMB
    fprintf('Running the comb - iteration #1\n')
    [~,~,~,pairwise,fmids] = combspectra(all_dat0_resamp,fcross,resamprate,parms,ifplot,mtmspecss);

    %% ------------ AMP/PHI QUALITY CONTROL METRICS --------------------
    % Calculate certain metrics to help work out if results can be trusted
    % This includes varying the amp2phiwt to see consistency between
    % phase/amplitude spectra
    fprintf('Testing for focusing\n')
    [gdofgd,focusscore] = amp2phi_QC_focusing(parms,pairwise,fmids,ifplot,{eqar(indgd).sta});

    if sum(gdofgd)<3 % <3 good events at all. save eqar wiped and move on
        fprintf('****************\n****************\n <3 good events without hint of focusing\n')
        if ifsave
                fprintf('Wiping all comb and saving eqar\n')
                save(arfile,'eqar')
        end
        fprintf('Moving on...\n')        
        continue
    end

    % subset to new "gd"
    indgd = indgd(gdofgd);
    fcross = fcross(gdofgd,:);
    all_dat0_resamp = all_dat0_resamp(:,gdofgd);
    mtmspecss.specs = mtmspecss.specs(:,gdofgd);

    %% ROUND 2 - RUN THROUGH COMB
    fprintf('Running the comb - iteration #2\n')
    [~,~,~,pairwise,fmids] = combspectra(all_dat0_resamp,fcross,resamprate,parms,ifplot,mtmspecss);


    %% pull out indiv
    Amat = pairwise.As;
    phimat = pairwise.phis;
    wtmat = double(pairwise.inds).*pairwise.wts;

    
    %% SAVE PAIRWISE SPECTRA 
    fprintf('%.0f Amp measurements\n',numel(Amat));
    sts = {datinfo(indgd).sta};
    if ifsave
        save(sprintf('RESULTS/PAIRSPECS/%.0f_pairspecs_comb_%s%s',evinfo.orids(ie),phase,component),'pairwise','sts','fmids','parms')
    end

    %% All in one inversion (if warranted)
    if all(all(wtmat==0)) % don't invert (and don't save results) if all junk
        delta_tstar_pref = nan(length(indgd),1);
        delta_tstar_std_pref = nan(length(indgd),1);
        delta_T_pref = nan(length(indgd),1);
        A0_pref = nan(length(indgd),1);
        alpha_pref = nan;
    else
    
        if isempty(parms.inv.alpha)
            %% ------------------  TEST FOR ALPHAS  ------------------
            test_alphas = [0:0.05:0.9];
    
            [ delta_tstar_pref,delta_T_pref,A0_pref,alpha_pref,alpha_VR,~,~,~,delta_tstar_std_pref ] ...
                = calc_fdependent( Amat,phimat,fmids,test_alphas,wtmat,parms.inv.amp2phiwt,parms.inv.opt,['Orid ',num2str(evinfo.orids(ie))] );
            close all
    
        else
            %% ------------------ DO ALL-IN-ONE INVERSION  ------------------
            [ delta_tstar_a0,delta_T_a0,A0_a0,~,~,~,~,~, delta_tstar_std_a0] ...
            = calc_fdependent( Amat,phimat,fmids,parms.inv.alpha,wtmat,parms.inv.amp2phiwt,parms.inv.opt,['Orid ',num2str(evinfo.orids(ie))],ifplot );
        
            % Assign a0 values
            delta_tstar_pref = delta_tstar_a0;
            delta_tstar_std_pref = delta_tstar_std_a0;
            delta_T_pref = delta_T_a0;
            A0_pref = A0_a0;
            alpha_pref = parms.inv.alpha;
            % parms.inv.alpha = 0;
        end
    end
    indnan = indgd(isnan(delta_tstar_pref));


    %% ---------------------- STORE RESULTS -----------------------
    % STORE RESULTS
    fprintf('Recording results in arrival structure...')
    % prep eqar to receive new fields
    eqar(1).dtstar_comb = []; 
    eqar(1).dtstarstd_comb = []; 
    eqar(1).dT_comb = []; 
    eqar(1).A0_comb = []; 
    eqar(1).alpha_comb = []; 
%     eqar(1).stds_comb = []; 
    eqar(1).par_dtstar_comb = parms;
    eqar(1).focusscore = [];
                      
    eqar(indgd) =  dealto(eqar(indgd),'dtstar_comb',delta_tstar_pref);
    eqar(indgd) =  dealto(eqar(indgd),'dT_comb',delta_T_pref);
    eqar(indgd) =  dealto(eqar(indgd),'A0_comb',A0_pref);
    eqar(indgd) =  dealto(eqar(indgd),'dtstarstd_comb',delta_tstar_std_pref);
    eqar(indgd) =  dealto(eqar(indgd),'alpha_comb',alpha_pref);
%     eqar(indgd) =  dealto(eqar(indgd),'stds_comb',std_dtstar_comb);
    eqar(indgd) =  dealto(eqar(indgd),'par_dtstar_comb',parms);
    eqar(indgd) =  dealto(eqar(indgd),'focusscore',focusscore(gdofgd));

	indbd = setdiff(1:length(eqar),indgd);
    if ~isempty(eqar(indbd))
        eqar(indbd) =  dealto(eqar(indbd),'dtstar_comb',nan);
        eqar(indbd) =  dealto(eqar(indbd),'dT_comb',nan);
        eqar(indbd) =  dealto(eqar(indbd),'A0_comb',nan);
    end

    %% while testing....
%     plot_ATTEN_TandF_domain_COMB_EAR( eqar(indgd) )
%     test_amp2phi_QC_focusing
%     return    
    
    
	%% -------------------------- PLOTS ---------------------------
    if ifplot
%       compare_dtstar_absAmp
        if strcmp(dbname,'EARdb')
            plot_ATTEN_TandF_domain_COMB_EAR( eqar(indgd) )
        else
            plot_ATTEN_TandF_domain_COMB( eqar(indgd) )
        end
    end % ifplot

%% -------------------------- SAVE ---------------------------
    if ifsave
    % SAVE
    save(arfile,'eqar')
    fprintf(' saved\n')
    % RECORD DTSTARCALC IN DATINFO 
    [datinfo.comb] = deal(false);
	[datinfo(indgd).comb] = deal(true); 
    [datinfo(indnan).comb] = deal(nan);
    save(datinfofile,'datinfo')
    
    fprintf(' STA  CHAN  NEZ  resp  tilt  comp  xcor  specR  comb\n')
	for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f     %1.0f     %1.0f     %1.0f     %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).rmtilt,datinfo(is).rmcomp,datinfo(is).xcor,datinfo(is).specR,datinfo(is).comb); end
    
    end

    
    toc  
end % loop on orids

cd '/Users/zeilon/Dropbox/MATLAB/AttenBody';
% results_PARSE