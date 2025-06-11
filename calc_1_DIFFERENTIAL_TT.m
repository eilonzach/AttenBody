% Script to cycle through events and calculate differential travel times
% for all stations. 
% 
% This script will plot a record section that allows you to select a window
% of data to cross-correlate and calculate differential travel times. You
% are then asked to select a nominal phase onset time for the corrected
% arrivals. 
% 
% Z. Eilon 2016
% 
% Commands:
%       x   select window (beginning and end)
%       f	change filter
%       s   toggle display of snr
%       y	toggle between gcarc distance and simply station order
%       d	delete station with trace beneath cursor
%       u	flip sign of station with trace beneath cursor
%       c	put on cheat-sheet of phases
%    left   prev event
%   right 	next event 
%      up 	increase amplitudes
%    down 	decrease amplitudes
%     ESC   stop

clear

% % project details
dbname = 'EARdb';
dbdir = '/Users/zeilon/Dropbox/Work/EARdb/'; % include final slash

% project details
% dbname = 'FRES_PILOT';
% dbdir = '~/Dropbox/Work/FRES_PILOT/'; % include final slash


%% parameters
phase = 'SKS'; % available are 'P,PP,S,PKP,PKS,SKS';
component = 'R'; %'Z', 'R', or 'T'
resamprate = 10 ; % new, common sample rate
% Do filtfs [12 1] for S, and [5 1] for P
% filtfs = 1./[2 0.25]; % P [flo fhi] = 1./[Tmax Tmin] in sec 
filtfs = 1./[12 1]; % S [flo fhi] = 1./[Tmax Tmin] in sec 
taperx = 0.2;
npoles = 2;
datwind = [-200 200]; % window of data in eqar structure
prelimwind = [-50 40]; % usually [-40 40], but may need to expand for longer periods
conflimwind = [-20 30]; % usually [-10 3y0], but may need to expand for longer periods
% prelimwind = [-25 15];
xcorlagmax = 20;
gcdistmax = 12;
mingcarc = 30;
maxgcarc = 130;
acormin = 0.75;
overwrite = false;
minNtrace = 4;

% parameters for the multitaper spectrogram of the phases
WinSize = 20;
FreqResolution=1/10;
StepSize= 0.25;   
% parameters for the multi-freq analysis of the phases
Tlos = [12 9   8  8  8  7 2 1 2 1.6666666 1 ];
This = [50 45 40 30 25 20 5 5 4 3.3333333 2]; % must be same length as above
% phase names for the cheat sheet (predicted arrival times
cheatsheetphases = {'P','S','PP','SKS','PKS','SKKS'};

%% Progress notes
% FRESPILOT -  S(T) got to 1863. P(Z) got to 1920. SKS(R) got to 1847

% EARdb_newpass_2024-02-14 finished P(Z). 
% EARdb_newpass_2024-02-21 finished PP(Z). 
% EARdb_newpass_2024-02-28 finished S(T). 
% EARdb_newpass_2024-02-28 finished PKS(R). 
% EARdb_newpass_2024-02-28 finished SKS(Z). 
% grabbed for most stations

%% EVT limits
firstev = 1; 
% time bounds for events to consider - ignore before startdate, or after enddate
startdate = '2014-03-01'; % format 'YYYY-MM-DD' % FOR EARdb new stations
% startdate = '1000-01-01'; % format 'YYYY-MM-DD' % for all
enddate   = '2025-01-01'; % format 'YYYY-MM-DD'
% magnitude bounds for events to consider - INCLUDE the limiting values
minmag = 6.2; % in Mw, set to 5 for all
maxmag = 7.2; % in Mw, set to 9 for all
% NOTE!! Also distance bounds for every phase, defined in the 
% "GRAB DATA IN RIGHT FORMAT" section below

% if you don't like any stations in particular, list them here to ignore
% their data in the record section and xcorr
% manualkillstas = {'C05E'};
manualkillstas = {''};



%% Preliminaries
wd = pwd;
addpath('matguts')
cd(dbdir);
run([dbdir,dbname,'_startup.m']);
run('~/Dropbox/MATLAB/startup.m')
spd = 24*60*60; % seconds per day - to convert time in seconds to serial time

%% Command list for TT picker
% 'x' = select each side of xcor window
% 'y' = change the y order (by distance or by simply in order)
% 'f' = change the filter
% 'd' = manually delete a station
%% 

%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

% event details
load([infodir,'/events'],'evinfo'); 
% station details
load([infodir,'/stations'],'stainfo');

% some default parameters
SNRdisp = false;
phasedisp = false;

ie = firstev; % up to 364
ied = 1; %incrementer - leave as 1
while 1 %evinfo.norids % loop on orids 
    tic
    close all, clear('x1','x2');
    if ie > evinfo.norids, error('ALL DONE!\n'); end
    orid = evinfo.orids(ie);
    fprintf('\n Orid %.0f %s \n\n',orid,evinfo.evtimes_IRISstr{ie}(1:end-4))
    % name files and directories
    evdir = [num2str(orid,'%03d'),'_',evinfo.datestamp(ie,:),'/'];
    datinfofile = [data_eqar_dir,evdir,'_datinfo_',phase];
    arfile      = [data_eqar_dir,evdir,'_EQAR_',phase,'_',component];

    % ignore outside date or magnitude bounds
    if  evinfo.evtimes(ie) < datenum(startdate) || evinfo.evtimes(ie) > datenum(enddate) ...
     || evinfo.evmags(ie) < minmag || evinfo.evmags(ie) > maxmag
        ie = ie + ied;
        nextev = 1;        
        continue; 
    end
   
    % check files exist
    if ~exist([datinfofile,'.mat'],'file'), fprintf('No station mat files for this event\n'); ie = ie+ied; continue, end
    if ~exist([arfile,'.mat'],'file'), fprintf('No arfile for this event + phase\n'); ie = ie + ied; continue, end
    
    % load files
    load(datinfofile) % loads datinfo stucture
    load(arfile)      % loads eqar structure
    
    % options to skip
    if isempty(datinfo), fprintf('No station mat files for this event\n'); ie = ie + ied; continue, end
    if isfield(datinfo,'xcor')
        if any([datinfo.xcor]==true)
            if ~overwrite
                yn = input('Xcor already done - overwrite? [y/n] ','s'); 
                if ~strcmp(yn,'y'), fprintf('skipping...\n'), ie = ie + ied; continue, end
            end
        end
    end
    if all([eqar.pred_arrT]==0), fprintf('No %s arrivals for this event\n',phase), ie = ie + ied; continue, end
    
    % account for different sample rate in eqar compared to resamprate
    if unique([eqar.samprate])~=resamprate
        fprintf('!! EQAR samprate is %.0f, whereas asking for resamprate = %.0f\n',unique([eqar.samprate]),resamprate)
        ynsr = input('Would you like to switch to EQAR samprate? ','s');
        if any(strcmp(ynsr,{'y','Y','yes'})) 
            fprintf('Okay, resamprate now %.0f\n',unique([eqar.samprate]))
            resamprate = unique([eqar.samprate]);
        else
            ie = ie + ied; continue, 
        end
    end

%     % skip if outside appropriate distance range for the phase
%     if strcmp(phase,'P') && mean([eqar.gcarc])>92, ie = ie + ied; continue, end
%     if strcmp(phase,'S') && mean([eqar.gcarc])>95, ie = ie + ied; continue, end
%     if strcmp(phase,'SKS') && mean([eqar.gcarc])<85, ie = ie + ied; continue, end
% 	if strcmp(phase,'SKP') && mean([eqar.gcarc])<129, ie = ie + ied; continue, end

    %% ------------------ GRAB DATA IN RIGHT FORMAT ------------------
    % prep data structures
    nstas = length(datinfo);
    wlen = diff(prelimwind)*resamprate*(1+4*taperx); % window length, accounting for taper+padding
    all_dat0  = zeros(resamprate*diff(datwind),nstas);
    dT_exist = nan(nstas,1);
    
    % LOOP ON STAS
    for is = 1:nstas
        fprintf('Station %s... ',datinfo(is).sta)        % APPLY DATA QUALITY CONDITIONS
        if isempty(eqar(is).(['dat' component])),fprintf('no data for this component\n'); continue; end   
        
        % skip if outside appropriate distance range for the phase        
        switch phase
            case 'P',   if eqar(is).gcarc > 99  || eqar(is).gcarc<30, fprintf('bad distance\n'),continue; end
            case 'PP',  if eqar(is).gcarc > 170 || eqar(is).gcarc<60, fprintf('bad distance\n'),continue; end
            case 'S',   if eqar(is).gcarc > 84  || eqar(is).gcarc<30, fprintf('bad distance\n'),continue; end
            case 'SKS', if eqar(is).gcarc > 132 || eqar(is).gcarc<88, fprintf('bad distance\n'),continue; end
            case 'SKKS',if eqar(is).gcarc > 132 || eqar(is).gcarc<91, fprintf('bad distance\n'),continue; end
            case 'PKS', if eqar(is).gcarc > 140 || eqar(is).gcarc<128, fprintf('bad distance\n'),continue; end
        end
        
        % GRAB DESIRED COMPONENT
        all_dat0(:,is) = eqar(is).(['dat',component])(:);

        % GRAB existing pick, if available
        if isfield(eqar,'abs_arrT')
            dT_exist(is) = eqar(is).abs_arrT;
        else
            dT_exist(is) = nan;
        end
        
        fprintf('got data\n')
    end % loop on stas
    
    % CLEAN DATA
    cp = struct('samprate',resamprate,'pretime',-datwind(1),...
            'prex',-prelimwind(1),'postx',prelimwind(2),...
            'taperx',taperx,'fhi',filtfs(2),'flo',filtfs(1),'npoles',npoles,'norm',1);
    [ all_datwf,all_datf,all_datwc,all_datc,~,ttws,tts ] = data_clean( all_dat0,cp );

    % MARK ZERO DATA TRACES
    indgd = 1:size(eqar); % << important variable - indgd are all the stations still in play...
    indgd(mean(abs(all_dat0(:,indgd)))<1e-16)     = []; % kill zero traces
    indgd(mean(abs(all_datwf(:,indgd)))==0)     = []; % kill zero traces
    indgd(isnan(mean(abs(all_datwf(:,indgd))))) = []; % kill nan traces
    for ik = 1:length(manualkillstas),indgd(strcmp({eqar(indgd).sta},manualkillstas{ik})) = []; end % kill bad stas traces
    
    %I think redundant after the above.... kill certain phases if likely to overlap with others
    if strcmp(phase,'SKKS') ,indgd([eqar(indgd).gcarc]<90) = []; end % kill close traces
    if strcmp(phase,'SKS') ,indgd([eqar(indgd).gcarc]<88) = []; end % kill close traces
    if strcmp(phase,'S') ,indgd([eqar(indgd).gcarc]>84) = []; end % kill far S traces

    % Decide on whether EQ is worth pursuing, with enough good traces
    if length(indgd) < minNtrace, fprintf('NOT ENOUGH GOOD TRACES/ARRIVALS, skip...\n'); ie = ie + ied; continue, end
    
    fprintf('Orid %.0f %s \n',orid,evinfo.evtimes_IRISstr{ie}(1:end-4))
	
    %% tauptime to help picking
    phasedb = mkcheatsheet(cheatsheetphases,[eqar(indgd).gcarc],evinfo.edeps(ie));
    % take off times of current phase to make differential times
    phasedb.dphasetimes = phasedb.phasetimes - phasedb.phasetimes(:,strcmp(phasedb.phases,phase))*ones(1,phasedb.Nph);
  
%% ------------------ PICK XCOR WINDOW, RECLEAN ------------------
    figxc = figure(53); clf(figxc); set(figxc,'position',[710 1 1626 1344])
    
%% ------------------ PICK TRAVEL TIMES ------------------

    % some default parameters
    cont=0;
    nextev=0;
    ydopt = 1;
    ampf = 1;
    autosnwin = 1;

    while nextev==0
    while cont==0

        %% Clear and set up the axes for the plotting
        try delete(ax1); end
        set(figxc,'position',[710 1 1626 1344])

        normf = ampf*max(max(abs(all_datwf(:,indgd))));
        % y option - distance or order
        if ydopt==1
            yds = [eqar(indgd).gcarc]; 
        else
            yds = 1:length(indgd);
        end
        yext = [min(yds),max(yds)];
        
        %make axes for plotting
        ax1 = axes(figxc,'pos',[0.05 0.07 0.3 0.85]); hold on
        cols = {'b','k'};

        % plot on the traces
        clear htr
        for ig = 1:length(indgd)
            is = indgd(ig);
            colr = cols{isodd(ig)+1};
            htr(ig) = plot(ax1,ttws,all_datwf(:,is)/normf/(10/diff(yext)) + yds(ig),'-','LineWidth',0.8,'color',colr);
            text(ax1,prelimwind(2)+1,yds(ig),datinfo(is).sta,'fontweight','bold','fontsize',20)

            % add previous picks if available
            plot(ax1,(dT_exist(is) - eqar(is).pred_arrT)*spd*[1;1],(diff(yext)/50)*[-1,1] + yds(ig),'-','LineWidth',1.5,'color',[1 0 0 0.5])
        end

        title(ax1,sprintf('PICK XCOR WINDOW\nORID %u, seaz %.0f, gcarc %.0f\ndep %.0fkm, Ms %.1f',...
            orid,mean([eqar.seaz]),mean([eqar.gcarc]),evinfo.edeps(ie),evinfo.evmags(ie)),'Fontsize',18)
        axis(ax1,[prelimwind min(yds)-0.06*diff(yext) max(yds)+0.06*diff(yext) ])
        set(ax1,'box','on','fontsize',16);
        try
            hp(1) = plot(ax1,[x1,x1],axlim(ax1,[3,4]),'--r','Linewidth',1.5);
            hp(2) = plot(ax1,[x2,x2],axlim(ax1,[3,4]),'--r','Linewidth',1.5);
        end

        % calculate some things re the filter
        centfreq = (filtfs(1)+filtfs(2))/2;
        centperiod=1/centfreq;
        if autosnwin
            nwinlim = [-3*centperiod , -1*centperiod];
            swinlim = [-1*centperiod , 1.5*centperiod];
        end


        %% SNR calculation and plot
        % if that flag is on
        % windowing 2.5 periods around the signal
        if SNRdisp == 1
            SNR = nan(size(all_datf,2),1);
            noisewin = find(nwinlim(1) < tts & tts <  nwinlim(2)); 
            sigwin   = find(swinlim(1) < tts & tts <  swinlim(2)); 

            % loop through good stations, calc SNR
            for ig = 1:length(indgd)
                is = indgd(ig);
               
%                 SNR(is) = (rms(all_dastf(sigwin,is))/rms(all_datf(noisewin,is)))^2;  
                SNR(is) = abs(maxab(all_datf(sigwin,is)))*rms(all_datf(sigwin,is))...
                          /...
                          rms(all_datf(noisewin,is))^2;  
            
                % plot the value SNR 
                text(ax1,prelimwind(1)+2,yds(ig),['SNR: ' num2str(SNR(is))],...
                    'fontsize',16,'VerticalAlignment','middle')

                colr = colour_get(SNR(is),max(SNR),min(SNR),redblue);
            end
            % now colour the traces on the plot accordingly
            for ig = 1:length(indgd)
                is = indgd(ig);
                colr = colour_get(log(SNR(is)),log(max(SNR)),log(min(SNR)),flipud(0.9*redblue).^2.5);
%                 colr = colour_get(SNR(is),max(SNR),min(SNR),flipud(0.9*jet).^3);
                set(htr(ig),'color',colr)
            end
            % plot histgram of the SNRs
            figsnr = figure(98);clf(figsnr),set(figsnr,'pos',[42 867 614 405])
            SNRbins = [0,0.5,0.8,1,1.5,2,3,4,5,10,20,30,40,50,75,100,150,200];
            hist(SNR,SNRbins),
            set(gca,'xscale','log','xlim',[0 200],'linewidth',1.5,'box','on','FontSize',15,'xtick',SNRbins)
            xlabel('SNR','FontSize',17,'FontWeight','bold')
            ylabel('counts','FontSize',17,'FontWeight','bold')

            % back to reality.
            figure(figxc);
        end

        %% add cheat sheet phases, if warranted
        if phasedisp==1 && ydopt==1
            % plot phase times, having subtracted off target time
            plot(ax1,phasedb.dphasetimes,phasedb.gcdists,':r','linewidth',1.2);
            % name the phases
            for ip = 1:phasedb.Nph
                gcplt = axlim(ax1,[3 4])*[0.98 0.02]'; % 1/25 of way up the plot
                ttplt = interp1(phasedb.gcdists,phasedb.dphasetimes(:,ip),gcplt);
                text(ax1,ttplt+1,gcplt,phasedb.phases{ip},'color','r','fontsize',15)
            end
        end
        



        %% Option to edit filter
        % work on figure 
        axes(ax1);
        [x,y,b] = ginput(1);
        if isempty(b), cont = 1; continue; end

        %% ================ PROCESS THE KEYBOARD INPUT ================
        switch b 

            case 'x'  % select window
                try delete(hp); end
                x1=x;
                hp(1) = plot(ax1,[x1,x1],axlim(ax1,[3,4]),'--r','Linewidth',1.5);
                [x2,y2] = ginput(1);
                hp(2) = plot(ax1,[x2,x2],axlim(ax1,[3,4]),'--r','Linewidth',1.5);
                xcorwind = sort(round_level([x1;x2],1./resamprate));

                % curtail filter if needed
                if diff(xcorwind) < 1/filtfs(1)
                    uiwait(warndlg('Window  length  shorter  than longest  period  in  filter... curtailing longest period','Window nüdge','modal'))
                    filtfs(1) = 1./diff(xcorwind);
                end

                % CLEAN DATA
                cp = struct('samprate',resamprate,'pretime',-datwind(1),...
                        'prex',-xcorwind(1),'postx',xcorwind(2),...
                        'taperx',taperx,'fhi',filtfs(2),'flo',filtfs(1),'npoles',npoles,'norm',1);
                [ xcor_datwf,~,xcor_datc,~,~,xcor_ttws,xcor_ttsc ] = data_clean( all_dat0,cp );

                % MARK ZERO DATA TRACES
                indgd(mean(abs(all_dat0(:,indgd)))<1e-16)     = []; % kill zero traces
                indgd(mean(abs(xcor_datwf(:,indgd)))==0)     = []; % kill zero traces
                indgd(isnan(mean(abs(xcor_datwf(:,indgd))))) = []; % kill nan traces
                if length(indgd) < 2
                    fprintf('NO GOOD TRACES/ARRIVALS, skip...\n')
                    [datinfo.xcor] = deal(false);
                    save(datinfofile,'datinfo');
                    ie = ie + ied;
                    cont=1;
                    continue
                end

                % Plot windowed and treated data
                try cla(ax2); end
                ax2 = axes(figxc,'pos',[0.38 0.07 0.27 0.85]); cla(ax2), hold on
                for ig=1:length(indgd)
                    is = indgd(ig);
%                     plot(ax2,xcor_ttws, xcor_datc(:,is)./max(abs(xcor_datc(:,is)))+ig,'k','Linewidth',1);
                    plot(ax2,xcor_ttws, xcor_datwf(:,is)./max(abs(xcor_datwf(:,is)))+ig,'r','Linewidth',1.5);
                    axis(ax2,[xcor_ttws(1) xcor_ttws(end) 0 ig+1]); 
                    text(ax2,xcor_ttws(end)+1,ig,datinfo(is).sta,'fontweight','bold','fontsize',20)
                    title(ax2,'Cleaned, filtered, windowed','Fontsize',18)
                end % loop on stas
                set(ax2,'box','on','ytick',[],'fontsize',16);


            case 'f' % change filter

                flt = inputdlg({'Enter  max  period','Enter  min  period','Enter  N  filter  poles'},...
                               'Filter details',1,{num2str(1/filtfs(1)),num2str(1/filtfs(2)),num2str(npoles)});
                flt = str2num(char(flt));
                filtfs(1) = 1/flt(1); % set high-pass (low-f corner)
                filtfs(2) = 1/flt(2); % set low-pass (high-f corner)
                npoles = flt(3);      % set N poles

                % re-filter data
                cp = struct('samprate',resamprate,'pretime',-datwind(1),...
                        'prex',-prelimwind(1),'postx',prelimwind(2),...
                        'taperx',taperx,'fhi',filtfs(2),'flo',filtfs(1),'npoles',npoles,'norm',1);
                [ all_datwf,all_datf,all_datwc,all_datc,~,ttws,tts ] = data_clean( all_dat0,cp );

                indgd(mean(abs(all_dat0(:,indgd)))<1e-16)     = []; % kill zero traces
                indgd(mean(abs(all_datwf(:,indgd)))==0)     = []; % kill zero traces
                indgd(isnan(mean(abs(all_datwf(:,indgd))))) = []; % kill nan traces

            case 'y'  % toggle between gcarc distance and simply station order

                ydopt = mod(ydopt+1,2);

            case 's'  % toggle on/off display of estimated SNR
                
                SNRdisp = ~SNRdisp;

            case 'v' % make stacked spectogram + other diagnostics  
                
                Make_Diagnostic_Figure_ze    
                figure(figxc);

            case 'c'  % toggle on/off display of cheat-sheet phases
                
                phasedisp = ~phasedisp;

            case 'q' % quality control by SNR - choose SNR cutoff

                snrqc = inputdlg({'Enter minimum SNR to keep'},'SNR-to-keep',1);
                if isempty(snrqc{:}), continue; end
                indgd( SNR(indgd) < eval(snrqc{:}) ) = []; % kill SNRs lower than cutoff

            case 'w' % choose signal and noise windows
                
                if autosnwin % if was automatic, now manual
                    title('FIRST CHOOSE SIGNAL WINDOW')
                    % choose signal window
                    [xs1,~] = ginput(1);
                    hsw(1) = plot(ax1,[xs1,xs1],axlim(ax1,[3,4]),'--b','Linewidth',1.5);
                    [xs2,~] = ginput(1);
                    hsw(2) = plot(ax1,[xs2,xs2],axlim(ax1,[3,4]),'--b','Linewidth',1.5);
                    swinlim = sort(round_level([xs1;xs2],1./resamprate));
                    title('NOW CHOOSE START OF NOISE WINDOW')
                     % choose noise window
                    [xn1,~] = ginput(1);
                    hnw(1) = plot(ax1,[xn1,xn1],axlim(ax1,[3,4]),'--r','Linewidth',1.5);
                    nwinlim = sort(round_level([xn1;swinlim(1)],1./resamprate));
                end
                % flip from auto to manual
                autosnwin = ~autosnwin;

                
            case 'd'  % delete station with trace beneath cursor

                indgd(mindex(yds,y)) = [];

            case 'u'  % flip sign of station with trace beneath cursor

                all_datwf(:,indgd(mindex(yds,y))) = -all_datwf(:,indgd(mindex(yds,y)));
                all_dat0(:,indgd(mindex(yds,y))) = -all_dat0(:,indgd(mindex(yds,y)));
                eqar(indgd(mindex(yds,y))).(['dat',component]) = -eqar(indgd(mindex(yds,y))).(['dat',component])(:);

            case 28 % left arrow
                break

            case 29  % right arrow
                break

            case 31 % up arrow
                ampf = ampf*1.4; % increase amplitudes

            case 30  % down arrow
                ampf = ampf/1.4; % decrease amplitudes

            case 27  % ESC key
                break

        end


        % move on to xcorr if ENTER selected
        if isempty(b), cont = 1; end
    
    end % continue to process 
    
    % Salsipuede
	if b==27  % ESCape if escaped
        fprintf('Skipping event on user request...\n')
        [datinfo.xcor] = deal(false);
        save(datinfofile,'datinfo');
        ie = ie + ied; % go forward one event
        nextev = 1;
        continue
    elseif b==28 % last event 
        fprintf('Previous event...\n')
        ied = -1;
        ie = ie+ied; % go back one event
        nextev = 1;
        continue
    elseif b==29 % next event 
        fprintf('Next event...\n')
        ied = 1;
        ie = ie + ied; % go forward one event
        nextev = 1;
        continue
    end
    
    %% ------------------ ITERATIVELY XCOR ------------------
    % removing chans with low acor
    fprintf('Cross-correlating + iteratively removing low acor traces\n')
    acor = zeros(size(eqar));
    
    acor_all = zeros(size(indgd));
    
    iter = 1;
    while any(acor < acormin)
        % DO CROSS CORRELATION
%         [dcor,dcstd,dvcstd,acor]=xcortimes(xcor_datwf(:,indgd),1./resamprate,xcorwind(1),xcorlagmax,0);
        [dcor,dcstd,dvcstd,acor]=...
            xcortimes_ze(xcor_datwf(:,indgd),1./resamprate,'t0',xcorwind(1),'lagmax',xcorlagmax,'ifwt',1,...
            'gclim',gcdistmax,'latlon',[[eqar(indgd).slat]',[eqar(indgd).slon]']);
        
        % kill nan dcors - todo far away/no constraint
        indgd(isnan(dcor)) = [];
        dcstd(isnan(dcor)) = [];
        dvcstd(isnan(dcor)) = [];
        acor(isnan(dcor)) = [];
        dcor(isnan(dcor)) = [];
        
        acor_all(indgd) = acor;
        
        % DELETE LOW ACOR CHANS
        if iter == 1 
            indgd(acor < 0.1) = []; % first kill really bad traces
        elseif iter > 1
            indgd(acor<acormin) = []; % then kill low acor traces
        end           
        if length(indgd) < 2, break, end
        iter=iter+1;
    end
	if length(indgd) < 2
        fprintf('NO GOOD TRACES/ARRIVALS, skip...\n')
        [datinfo.xcor] = deal(false);
        save(datinfofile,'datinfo');
        nextev = 1;
        continue
    end
    
    %% ------------------ PLOTTING ------------------
    % RECORD SECTION TO CHECK
    try cla(ax3); end
    ax3 = axes(figxc,'pos',[0.68 0.07 0.27 0.85]); cla(ax3)
    set(ax3,'ytick',[],'box','on','fontsize',16,...
        'xlim',conflimwind); 
    hold on
    normf = max(max(abs(all_datwf(:,indgd))));
%     normfx = max(max(abs(xcor_datwf(:,indgd))));
    gcarcs = [eqar.gcarc];
    
    kk = 0;
    for ig = 1:length(indgd)
        is = indgd(ig);
        kk = kk+1;
%         plot(ttws                   ,all_datwf(:,is)/normf/2 + gcarcs(is),'k')
        plot(ax3,ttws-dcor(ig),all_datwf(:,is)/normf + kk,'b','LineWidth',1.5)
%         plot(ax3,xcor_ttws-dcor(ig),xcor_datwf(:,is)/normfx + kk,'g','LineWidth',1.5)
        text(ax3,axlim(ax3,2)+1,kk,datinfo(is).sta,'fontweight','bold')
    end

	% STACK raw traces to get overall arrtime
    stk = zeros(wlen+1,1);
    for ig = 1:length(indgd)
        is = indgd(ig);
        nshift = round(dcor(ig)*resamprate);
        % can only do this trick because padded with zeros

        datstk = [zeros(-2*nshift,1);all_datwf(abs(nshift)+1:end-abs(nshift),is);zeros(2*nshift,1)];
        stk = stk + datstk;
    end
    plot(ax3,ttws(1:wlen+1),stk./max(abs(stk)) + axlim(ax3,3),'r','LineWidth',1.5)
	text(ax3,ttws(end)+0.5,axlim(ax3,3),'*STACK*')
	xlabel(ax3,sprintf('Orid %.0f,  gcarc ~%.1f,  seaz ~%.1f',ie,mean([eqar.gcarc]),mean([eqar.seaz])),'Fontsize',22)
    ylim(ax3,[-1 1+length(indgd)])
   
    % PICK phase arrival
    conf = 0;
    while conf==0
        title('PICK PHASE ARRIVAL TIME (ESC or > TO SKIP)','FontSize',18)
        axes(ax3)
        [t0,~,b3] = ginput(1);
        switch b3
            case {1,'x'} % click or 'x'
            try delete(hl); end
            hold on 
            hl = plot([t0(1),t0(1)],[axlim(ax3,3),axlim(4)],'--k');
            title('ENTER TO CONFIRM','FontSize',18)
            axes(ax3)
            [~,~,b4] = ginput(1);
            if isempty(b4)
                conf=1;
            end
            
            case {27,29}  % ESC or forward key
                conf = nan;
                ied = 1;
                break
            case 28  % back key
                conf = nan;
                ied = 0;
                break
        end
    end
	if isnan(conf)
        [datinfo.xcor] = deal(false);
        save(datinfofile,'datinfo');
        ie = ie + ied;
        nextev = 1;        
        continue
    end

    
    %% ---------------------- STORE RESULTS -----------------------
    % STORE RESULTS
    fprintf('Recording results in arrival structure...')
    % prep eqar to receive new fields
    eqar(1).dT = []; eqar(1).acor = []; eqar(1).abs_arrT = []; eqar(1).par_dT = [];

    eqar(indgd) =  dealto(eqar(indgd),'dT',dcor);
    eqar(indgd) =  dealto(eqar(indgd),'acor',acor);
    eqar(indgd) =  dealto(eqar(indgd),'abs_arrT',[eqar(indgd).pred_arrT]' + dcor./spd + t0(1)./spd);
    eqar(indgd) =  dealto(eqar(indgd),'par_dT',cp);

	indbd = setdiff(1:length(eqar),indgd);
    if ~isempty(indbd)
        eqar(indbd) =  dealto(eqar(indbd),'dT',nan);
        eqar(indbd) =  dealto(eqar(indbd),'acor',nan);
        eqar(indbd) =  dealto(eqar(indbd),'abs_arrT',nan);
        eqar(indbd) =  dealto(eqar(indbd),'par_dT',nan);
    end

    
    pause(.1)


    % SAVE
    save(arfile,'eqar')
    fprintf(' saved\n')
    % RECORD XCOR IN DATINFO 
    [datinfo.xcor] = deal(false);
    [datinfo(indgd).xcor] = deal(true);
    save(datinfofile,'datinfo')
    
    fprintf(' STA  CHAN  NEZ  resp  tilt  comp  xcor\n')
	for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f     %1.0f     %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).rmtilt,datinfo(is).rmcomp,datinfo(is).xcor); end
    
    % move on to next event
    ied = 1;
    ie = ie + ied;
    nextev = 1;
    ied = 1;
    end 
    toc  
    
end % while processing