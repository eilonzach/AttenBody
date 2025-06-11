% cycle through events and calculate attenuation for all stations by
% running a comb of filters over each pair of stations' traces and
% calculating relative phase and differential amplitude at each frequency,
% before fitting the spectra simultaneously with curves whose slopes relate
% simply to delta-tstar
close all

% project details
% dbname = 'EARdb';
% dbdir = '/Users/zeilon/Dropbox/Work/EARdb/'; % include final slash


% dbname = 'AKTAdb';
% dbdir = '/Users/zeilon/Work/AKTA_Cristhian/AKTAdb/'; % include final slash

% % project details
dbname = 'EARdb';
dbdir = '/Users/zeilon/Dropbox/Work/EARdb/'; % include final slash

% project details
% dbname = 'FRES_PILOT';
% dbdir = '~/Dropbox/Work/FRES_PILOT/'; % include final slash


%% parameters
phase = 'P';
component = ['Z']; %'Z', 'R', or 'T'

overwrite  = false; % should be "true" if ifplotonly 
ifsave     = false;
ifplot     = 0; % 2 for the individual combs...
ifplotonly = false;

% datproc parms
arrT_pred_or_abs = 'abs'; % option to align traces by predicted (tauP) time or the absolute measured diffTT time
filtfs = 1./[40 0.5]; % [flo fhi] = 1./[Tmax Tmin] in sec
taperx = 0.1;
datwind = [-160 160]; % window of data in eqar structure
% note, inside fAxPspectra will pad this so as to not lose data due to taper
specwind = [-10 10]; % [-prex postx] both values are relative to arrival (prex is thus +ive)
snrmin = 4;

%% EVT limits
firstev = 1291; 
% time bounds for events to consider - ignore before startdate, or after enddate
startdate = '2000-03-01'; % format 'YYYY-MM-DD' % FOR EARdb new stations
% startdate = '1000-01-01'; % format 'YYYY-MM-DD' % for all
enddate   = '2025-01-01'; % format 'YYYY-MM-DD'
% magnitude bounds for events to consider - INCLUDE the limiting values
minmag = 6.2; % in Mw, set to 5 for all
maxmag = 7.2; % in Mw, set to 9 for all
% NOTE!! Also distance bounds for every phase, defined in the 
% "GRAB DATA IN RIGHT FORMAT" section below

%% fAxP parms
parms.fAxP.ifunwrap = 0;

parms.wind.pretime = -datwind(1); % seconds before the predicted arrival that the data time series starts. Should be a positive number
parms.wind.prex = -specwind(1); % seconds before the predicted arrival to take window. So if you want to start before the phase arrive, this should be positive. 
parms.wind.postx = specwind(2); % seconds after the predicted arrival to take window. So if you want to start after the phase arrive, this should be positive. 
parms.wind.taperx = taperx;
parms.wind.flo = filtfs(1);
parms.wind.fhi = filtfs(2);

parms.qc.minacor = 0.6;
parms.qc.maxphi = 5;
parms.qc.mtmRcomp = true; % option to compare comb and mtm estimates of As
parms.qc.phi2deriv = true; % option to weight by second derivative of phi not being huge
parms.qc.focus_threshold = 3; % threshold value for focus metric to flag up as bad
parms.qc.maxdist = 5; % maximum distance, in degrees, between stations to consider

parms.inv.amp2phiwt = 2; % ratio of upweight of amplitude vs. phase spectra fit

parms.inv.fmin_cskip = 0.1; 
parms.inv.fmin = 0.045; % caps maximum value otherwise set by fcrosslo
parms.inv.fmax = 1;    % caps minimum value otherwise set by fcrosshi
parms.inv.corr_c_skip = true; % correct for cycle-skip. 
parms.inv.ifwt = true; % weight data by xcorr value in fit of narrow-band filtered data between two stations
parms.inv.R2default = 0.9; % if zero, do not weight by R2 of each pairwise fit, Else, weight by multiple of this default.
parms.inv.opt = 1; %USE 1!   for calc_fdependent...  1 is all in one method; 2 is one-by-one method

parms.inv.alpha = 0.27; % assumed frequency dependency of attenuation - set to empty to test for this
% parms.inv.alpha = []; % assumed frequency dependency of attenuation - set to empty to test for this

if strcmp(dbname,'AKTAdb')
    datwind = [-180 90]; % window of data in eqar structure
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


for ie = firstev:evinfo.norids %evinfo.norids % 335:norids % loop on orids
%     if  mags(ie)<6.9, continue, end
    tic
    orid = evinfo.orids(ie);
    fprintf('\n Orid %.0f %s %s-%s\n\n',orid,evinfo.evtimes_IRISstr{ie}(1:end-4),phase,component)
    % name files and directories
    evdir = [num2str(orid,'%03d'),'_',evinfo.datestamp(ie,:),'/'];
    datinfofile = [data_eqar_dir,evdir,'_datinfo_',phase];
    arfile      = [data_eqar_dir,evdir,'_EQAR_',phase,'_',component];

    % ignore outside date or magnitude bounds
    if  evinfo.evtimes(ie) < datenum(startdate) || evinfo.evtimes(ie) > datenum(enddate) ...
     || evinfo.evmags(ie) < minmag || evinfo.evmags(ie) > maxmag
        continue; 
    end

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
    if isfield(datinfo,'fAxP')
        if any([datinfo.specR]==true)
            if ~ifplotonly
            if ~overwrite 
                yn = input('delta-tstar_fAxP already done - overwrite? [y/n] ','s'); 
                if ~strcmp(yn,'y'), fprintf('skipping...\n'), continue, end
            end
            end
        end
    end
    
    if ifplotonly
        indgd = ~isnan([eqar.dtstar_fAxP]);
        if ~any(indgd), continue; end
%         if ~any([eqar(indgd).slat]>12.5 & [eqar(indgd).slon]<42  & [eqar(indgd).slon]>39), continue; end
        plot_ATTEN_TandF_domain_fAxP_EAR(eqar(indgd))
        return
        continue
    end

    % reset all fAxP-related fields in eqar if saving
    if ifsave
        for fnames = {'dtstar_fAxP','dT_fAxP','A0_cfAxP','alpha_fAxP'}
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
        if isnan(eqar(is).dT),fprintf('nan xcor for this component\n'),continue; end
        try
            if isnan(eqar(is).dtstar_specR),fprintf('no specR-dtstar for this component\n'),continue; end
        catch
            warning('eqar seems to call specR method "dtstar"')
            if isnan(eqar(is).dtstar),fprintf('no specR-dtstar for this component\n'),continue; end
        end

        % SHIFT TRACES USING Predicted ARRIVAL TIME 
        % OR... 
        % SHIFT TRACES USING Observed diffTT ARRIVAL TIME
        att = (eqar(is).tt-eqar(is).([arrT_pred_or_abs,'_arrT']))*spd + dcTshft; % shift to since predicted arrival 
        att = round(att,4);
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
    %  isob = zeros(size(eqar)); for is = 1:length(eqar), isob(is) = ~isempty(which_OBS(eqar(is).sta)); end
    
    % ONLY USE GOOD TRACES
    indgd = 1:length(eqar);
    indgd(mean(abs(all_dat0(:,indgd)))==0)     = []; % kill zero traces
    indgd(mean(abs(all_dat0(:,indgd)))<1e-12)     = []; % kill zero traces
    indgd(isnan(mean(abs(all_dat0(:,indgd))))) = []; % kill nan traces
    try
        indgd(isnan([eqar(indgd).dtstar_specR])) = []; % kill traces without specR dtstar 
    catch
        warning('eqar seems to call specR method "dtstar"')
        indgd(isnan([eqar(indgd).dtstar])) = []; % kill traces without specR dtstar 
    end
    if ~isfield(eqar,'snr_wf'), if overwrite,keyboard; end, continue; end % skip if snr not even calculated
    indgd([eqar(indgd).snr_wf]<snrmin)                  = []; % kill low snr traces
    indgd([eqar(indgd).acor] < median([eqar(indgd).acor]) - 2*std([eqar(indgd).acor])) = []; % kill traces with acor that (somehow) passed QC but is actually more than 2std lower than median

    if parms.qc.maxdist
    gdofgd = true(size(indgd));
    % get slatlon
        slatlon = nan(length(eqar),2);
        for is = 1:length(eqar)
            if ~isempty(eqar(is).slat), slatlon(is,:) = [eqar(is).slat,eqar(is).slon]; end
        end
        for is = indgd
            d2all = distance(slatlon(is,:),slatlon(indgd,:));
            d2all = d2all(indgd~=is);
            if sum(d2all<parms.qc.maxdist)<3 % if no more than 2 sta within maxdist...
                gdofgd(indgd==is) = 0; 
            end
        end
        indgd = indgd(gdofgd);
        clear gdofgd
    end



    % only keep good stas
    all_dat0_gd = all_dat0(:,indgd);

    % kick out weird gain stations (log of rms/median(rms) > 2)
    indgd(abs(log(rms(all_dat0_gd)./median(rms(all_dat0_gd))))>2) = [];
    all_dat0_gd = all_dat0(:,indgd);

    if length(indgd) < 2, fprintf('NO GOOD TRACES/ARRIVALS, skip...\n'), 
            if ifsave && overwrite
                fprintf(['***********************************\n' ...
                    'Wiping all fAxP and saving eqar\n' ...
                    '***********************************\n'])
                save(arfile,'eqar')
            end
        continue
    end

    % get spectral ratios for comparison
    mtmspecss = struct('frq',eqar(indgd(1)).frq,'specs',[eqar(indgd).specs]);

    % get slatlon
    slatlon = [[eqar(indgd).slat]',[eqar(indgd).slon]'];
    
    %% GET station-by-station acceptable f range from crossing freqs
    fcross = [[eqar(indgd).fcrosslo]',[eqar(indgd).fcrosshi]'];
    fcross(fcross(:,1) < parms.inv.fmin,1) = parms.inv.fmin;
    fcross(fcross(:,2) > parms.inv.fmax,2) = parms.inv.fmax;

    % get sample rate
    samprate = unique([eqar(indgd).samprate]);
    if length(samprate)>2, error('seem to have multiple sample rates; need to resample to common rate'); end
    
    %% resamp - speeds up - NOT NEEDED
%     fprintf('resampling to %.0f Hz\n',resamprate);
%     tt0 = att(ja)';
%     tt1 = [tt0(1):1/resamprate:tt0(end)]';
%     nsamps = length(tt1); nstas = size(all_dat0_gd,2);
%     all_dat0_resamp = zeros(nsamps,nstas);
%     for is = 1:nstas
%         all_dat0_resamp(:,is) = interp1(tt0,all_dat0_gd(:,is),tt1);
%     end

    %% DO 2X! RUN THROUGH fAxP

    for itx = 1:2
        
        fprintf('Running the fAxP calculation - iteration #%.0f\n',itx)
        [~,~,~,pairwise,fmids] = fAxPspectra(all_dat0_gd,fcross,samprate,parms,ifplot,mtmspecss,[[eqar(indgd).slat]',[eqar(indgd).slon]']);
    
        %% ------------ AMP/PHI QUALITY CONTROL METRICS --------------------
        % Calculate certain metrics to help work out if results can be trusted
        % This includes varying the amp2phiwt to see consistency between
        % phase/amplitude spectra
        fprintf('Testing for focusing\n')
        [gdofgd,focusscore] = amp2phi_QC_focusing(parms,pairwise,fmids,1,{eqar(indgd).sta});
        
        if sum(gdofgd)<3 % <3 good events at all. save eqar wiped and move on
            break
        end

        % subset to new "gd"
        indgd = indgd(gdofgd);
        fcross = fcross(gdofgd,:);
        all_dat0_gd = all_dat0_gd(:,gdofgd);
        mtmspecss.specs = mtmspecss.specs(:,gdofgd);

    end

    if sum(gdofgd)<3 % <3 good events at all. save eqar wiped and move on
        fprintf('****************\n****************\n <3 good events without hint of focusing\n')
        if ifsave
                fprintf('Wiping all comb and saving eqar\n')
                save(arfile,'eqar')
        end
        fprintf('Moving on...\n')        
        continue
    end

    % do one last time to get pairwise that only has gdofgd in
    fprintf('Running the fAxP calculation - iteration #%.0f\n',itx+1)
    [~,~,~,pairwise,fmids] = fAxPspectra(all_dat0_gd,fcross,samprate,parms,ifplot,mtmspecss,slatlon);

    %% pull out indiv
    Amat = pairwise.As;
    phimat = pairwise.phis;
    wtmat = double(pairwise.inds).*pairwise.wts;

    
    %% SAVE PAIRWISE SPECTRA 
    fprintf('%.0f Amp measurements\n',numel(Amat));
    sts = {datinfo(indgd).sta};
    if ifsave
        save(sprintf('RESULTS/PAIRSPECS/%.0f_pairspecs_fAxP_%s%s',evinfo.orids(ie),phase,component),'pairwise','sts','fmids','parms')
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
    eqar(1).dtstar_fAxP = []; 
    eqar(1).dtstarstd_fAxP = []; 
    eqar(1).dT_fAxP = []; 
    eqar(1).A0_fAxP = []; 
    eqar(1).alpha_fAxP = []; 
%     eqar(1).stds_fAxP = []; 
    eqar(1).par_dtstar_fAxP = parms;
    eqar(1).focusscore_fAxP = [];
                      
    eqar(indgd) =  dealto(eqar(indgd),'dtstar_fAxP',delta_tstar_pref);
    eqar(indgd) =  dealto(eqar(indgd),'dT_fAxP',delta_T_pref);
    eqar(indgd) =  dealto(eqar(indgd),'A0_fAxP',A0_pref);
    eqar(indgd) =  dealto(eqar(indgd),'dtstarstd_fAxP',delta_tstar_std_pref);
    eqar(indgd) =  dealto(eqar(indgd),'alpha_fAxP',alpha_pref);
    eqar(indgd) =  dealto(eqar(indgd),'par_dtstar_fAxP',parms);
    eqar(indgd) =  dealto(eqar(indgd),'focusscore_fAxP',focusscore(gdofgd));

	indbd = setdiff(1:length(eqar),indgd);
    if ~isempty(eqar(indbd))
        eqar(indbd) =  dealto(eqar(indbd),'dtstar_fAxP',nan);
        eqar(indbd) =  dealto(eqar(indbd),'dT_fAxP',nan);
        eqar(indbd) =  dealto(eqar(indbd),'A0_fAxP',nan);
    end

    %% while testing....
%     plot_ATTEN_TandF_domain_fAxP_EAR( eqar(indgd) )
%     test_amp2phi_QC_focusing
%     return    
    
    
	%% -------------------------- PLOTS ---------------------------
    if ifplot
%       compare_dtstar_absAmp
        if strcmp(dbname,'EARdb')
            plot_ATTEN_TandF_domain_fAxP_EAR( eqar(indgd) )
        elseif strcmp(dbname,'FRES_PILOT')
            plot_ATTEN_TandF_domain_fAxP_FRES_PILOT( eqar(indgd) )
        else
            plot_ATTEN_TandF_domain_fAxP( eqar(indgd) )
        end
    end % ifplot

%% -------------------------- SAVE ---------------------------
    if ifsave
    % SAVE
    save(arfile,'eqar')
    fprintf(' saved\n')
    % RECORD DTSTARCALC IN DATINFO 
    [datinfo.fAxP] = deal(false);
	[datinfo(indgd).fAxP] = deal(true); 
    [datinfo(indnan).fAxP] = deal(nan);
    save(datinfofile,'datinfo')
    
    try
    fprintf(' STA  CHAN  NEZ  resp  tilt  comp  xcor  specR  comb  fAxP\n')
	for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f     %1.0f     %1.0f     %1.0f     %1.0f     %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).rmtilt,datinfo(is).rmcomp,datinfo(is).xcor,datinfo(is).specR,datinfo(is).comb,datinfo(is).fAxP); end
    catch % don't even have a comb for this one!
	for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f     %1.0f     %1.0f     %1.0f     %1.0f     0     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).rmtilt,datinfo(is).rmcomp,datinfo(is).xcor,datinfo(is).specR,datinfo(is).fAxP); end
    end

    end % on ifsave

    
    toc  
end % loop on orids

cd '/Users/zeilon/Dropbox/MATLAB/AttenBody';
% results_PARSE