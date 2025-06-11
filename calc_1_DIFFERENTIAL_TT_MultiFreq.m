% Script to cycle through events and calculate differential travel times
% for all stations. 
% 
% This script will plot a record section that allows you to select a window
% of data to cross-correlate and calculate differential travel times. You
% are then asked to select a nominal phase onset time for the corrected
% arrivals. 
% 
% Z. Eilon 2016
% A. Hariharan 2024 added Multifrequency capability
% 
% Commands:
%       x   select window (beginning and end)
%       f	change filter
%       y	toggle between gcarc distance and simply station order
%       d	delete station with trace beneath cursor
%       u	flip sign of station with trace beneath cursor
%    left   prev event
%   right 	next event 
%      up 	increase amplitudes
%    down 	decrease amplitudes
%     ESC   stop

%% Command list for TT picker
% 'x' = select each side of xcor window
% 'y' = change the y order (by distance or by simply in order)
% 'f' = change the filter
% 'd' = manually delete a station
% 'v' = make a diagnostic spectogram figure
% 'k' = Turn cluster analysis on/off
% 'c' = Zoom in/out on the window for picking, but remember traces are
% already tapered
% 'z' = Zoom in/out of the window for assigning phase arrival time
% 't' = Display SNR of a given trace on
%
clear

%% Preliminaries: Make sure the lines below are calibrated for your file system

% project details
dbname = 'EARdb';
dbdir = '/Users/zeilon/Dropbox/Work/EARdb/'; % include final slash

% wd = pwd;
% addpath('/Volumes/data/Anant/OldOrca_Tomo/AttenBody-main1/matguts')
% addpath('/Volumes/data/Anant/OldOrca_Tomo/AttenBody-main1/seis_tools')
% addpath('/Volumes/data/Anant/OldOrca_Tomo/myFUNCTIONS-master')
% addpath('/Volumes/data/Anant/OldOrca_Tomo')
% addpath('/Volumes/data/Anant/OldOrca_Tomo/multitaper_toolbox-master/matlab')
% dbname = 'OldORCA_Try2';
% dbdir = '/Volumes/data/Anant/OldOrca_Tomo/OldORCA/'; % include final slash
% run('/Volumes/data/Anant/OldOrca_Tomo/OldORCA_Try2/OldORCA_Try2_startup.m');
% data_eqar_dir = '/Volumes/data/Anant/OldOrca_Tomo/OldORCA_Try2/EQAR_NONOISE/';
% addpath('/Volumes/data/Anant/OldOrca_Tomo/AttenBody-main1')
% cd(dbdir); 

%% parameters
phase = 'P'; %S';%'P';
component = 'Z'; %'T';%'Z'; %'Z', 'R', or 'T'
resamprate = 10 ; % new, common sample rate
overwrite = false;

% Do filtfs [12 1] for S, and [5 1] for P
% filtfs = 1./[2 0.25]; % P [flo fhi] = 1./[Tmax Tmin] in sec 
filtfs = 1./[2.5 0.5]; % S [flo fhi] = 1./[Tmax Tmin] in sec 
taperx = 0.2;
npoles = 2;
datwind = [-200 200]; % window of data in eqar structure
prelimwind = [-80 80];
% prelimwind = [-25 15];
xcorlagmax = 20;
gcdistmax = 12;
mingcarc = 30;
maxgcarc = 130;
acormin = 0.75;
minNtrace = 4;

% parameters for the multitaper spectrogram of the phases
WinSize = 20;
FreqResolution=1/10;
StepSize= 0.25;     

clusteroff=0;
OverlayPhase=0;
mint=0.5;
hit=10;
Nwds=10;
SNRdisp=0;

% FRESPILOT -  S(T) got to 1863. P(Z) got to 1920. SKS(R) got to 1847
% EARdb (I think outdated) SKS got to 1044; S got to 1287; % need to do 160-335 again on all
firstev = 1; 
% time bounds for events to consider - ignore before startdate, or after enddate
startdate = '2014-03-01'; % format 'YYYY-MM-DD' 
enddate   = '2025-02-01'; % format 'YYYY-MM-DD'


%% Make set of period windows for bandpass filter
Tlos = [12 9   8  8  8  7 2 1 2 1.6666666 1 ];
This = [50 45 40 30 25 20 5 5 4 3.3333333 2];
Tmids = (Tlos+This)./2;
flos = 1./This;
fhis = 1./Tlos;
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

 

%% 

%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%

% event details
load([infodir,'/events'],'evinfo'); 
% station details
load([infodir,'/stations'],'stainfo');

ie = 29; % up to 364
ied = 1; %incrementer - leave as 1
for ie= ie:ied:500 % loop on orids 
    FCounter=1;
%     if  mags(ie)<6.9, continue, end
    tic
    close all, clear('x1','x2');
    orid = evinfo.orids(ie);
    fprintf('\n Orid %.0f %s \n\n',orid,evinfo.evtimes_IRISstr{ie}(1:end-4))
    % name files and directories
    evdir = [num2str(orid,'%03d'),'_',evinfo.datestamp(ie,:),'/'];
    datinfofile = [data_eqar_dir,evdir,'_datinfo_',phase];
    arfile      = [data_eqar_dir,evdir,'_EQAR_',phase,'_',component];

     % ignore outside date bounds
    if evinfo.evtimes(ie) < datenum(startdate) || evinfo.evtimes(ie) > datenum(enddate), 
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
    datinfo.rmtilt
    % options to skip
    if isempty(datinfo), fprintf('No station mat files for this event\n'); ie = ie + ied; continue, end
    if isfield(datinfo,'xcor')
        if any([datinfo.xcor]==true)
            if ~overwrite
                yn = input('Xcor already done - overwrite? [y/n] ','s'); 
                if strcmp(yn,'y')
                    if isfield(datinfo,'xcor')
                    JUNK=rmfield(datinfo,'xcor');
                    end
                    if isfield(datinfo,'Freqlolist')
                    JUNK=rmfield(datinfo,'Freqlolist');
                    end
                    if isfield(datinfo,'Freqhilist')
                    JUNK=rmfield(datinfo,'Freqhilist');    
                    end
                end
                if ~strcmp(yn,'y'), fprintf('skipping...\n'), ie = ie + ied; continue, end
            end
        end
    end
    if all([eqar.pred_arrT]==0), fprintf('No %s arrivals for this event\n',phase), ie = ie + ied; continue, end
    
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
    
    % LOOP ON STAS
    for is = 1:nstas
        fprintf('Station %s... ',datinfo(is).sta)        % APPLY DATA QUALITY CONDITIONS
        if isempty(eqar(is).(['dat' component])),fprintf('no data for this component\n'); continue; end   
        
        % skip if outside appropriate distance range for the phase        
        switch phase
            case 'P',   if eqar(is).gcarc > 95 || eqar(is).gcarc<30, continue; end
            case 'S',   if eqar(is).gcarc > 84 || eqar(is).gcarc<30, continue; end
            case 'SKS', if eqar(is).gcarc > 132 || eqar(is).gcarc<86, continue; end
            case 'PKS',   if eqar(is).gcarc > 140 || eqar(is).gcarc<128, continue; end
        end
        
        % GRAB DESIRED COMPONENT
        all_dat0(:,is) = eqar(is).(['dat',component])(:);
        
        fprintf('got data\n')
    end % loop on stas
    
    % CLEAN DATA
    cp = struct('samprate',resamprate,'pretime',-datwind(1),...
            'prex',-prelimwind(1),'postx',prelimwind(2),...
            'taperx',taperx,'fhi',filtfs(2),'flo',filtfs(1),'npoles',npoles,'norm',1);
    [ all_datwf,all_datf,all_datwc,all_datc,~,ttws,tts ] = data_clean( all_dat0,cp );

    % MARK ZERO DATA TRACES
    indgd = 1:size(eqar); % << important variable - indgd are all the stations still in play...
    %indgd(mean(abs(all_dat0(:,indgd)))<1e-16)     = []; % kill zero traces
%     indgd(mean(abs(all_datwf(:,indgd)))==0)     = []; % kill zero traces
%     indgd(isnan(mean(abs(all_datwf(:,indgd))))) = []; % kill nan traces
    for ik = 1:length(manualkillstas),indgd(strcmp({eqar(indgd).sta},manualkillstas{ik})) = []; end % kill bad stas traces
    
    if strcmp(phase,'SKS') ,indgd([eqar(indgd).gcarc]<86) = []; end % kill close traces
    
    if length(indgd) < minNtrace, fprintf('NOT ENOUGH GOOD TRACES/ARRIVALS, skip...\n'); ie = ie + ied; continue, end
    
    fprintf('Orid %.0f %s \n',orid,evinfo.evtimes_IRISstr{ie}(1:end-4))
	%% tauptime to help picking
    % tauptime('d',mean([eqar.gcarc]),'z',evinfo.edeps(ie))
    % 
    % %% Compile arrival times for different phases
    % MegaPhaseList= [];
    % for temptempcount = 1:length(eqar)
    %         distlist(temptempcount) = eqar(temptempcount).gcarc;
    %         tmpphaselist=extractfield(tauptime('d',([eqar(temptempcount).gcarc]),'z',evinfo.edeps(ie)),'phase');
    %         MegaPhaseList=[MegaPhaseList tmpphaselist];  
    % end
    % PhasesforthisEQ = unique(MegaPhaseList);
    % 
    %    for temptempcount = 1:length(eqar)
    %        currttimerun=tauptime('d',([eqar(temptempcount).gcarc]),'z',evinfo.edeps(ie));
    %          tmpphaselist=extractfield(currttimerun,'phase');
    %         thisphase = find(strcmp(tmpphaselist,phase)==1 | strcmp(tmpphaselist,lower(phase))==1);
    %         thisphase=thisphase(1);
    %         for phasechecker=1:length(PhasesforthisEQ)
    %             currphase=PhasesforthisEQ{phasechecker};
    %             matchdx = find(strcmp(currphase,tmpphaselist)==1);
    %              if length(matchdx) == 0
    %                  ArrivalTimeList(phasechecker,temptempcount) = NaN;
    %              else
    % 
    %             matchdx=matchdx(1);
    % 
    % 
    %             ArrivalTimeList(phasechecker,temptempcount) = currttimerun(matchdx).time-currttimerun(thisphase).time;
    %              end
    %         end
    % 
    % 
    %    end

%% ------------------ PICK XCOR WINDOW, RECLEAN ------------------
    figxc = figure(53); clf(figxc); set(figxc,'position',[710 1 1626 1344])
    
    %% PICK TRAVEL TIMES
    cont=0;
    nextev=0;
    ydopt = 1;
    ampf = 1;
    while nextev==0
    while cont==0
        try clear SNR; end
            centfreq = (filtfs(1)+filtfs(2))/2;
            centperiod=1/centfreq;
        try delete(ax1); end
        sigwin = find(ttws > -1*centperiod & ttws < 1.5*centperiod)
        % cluster analysis as very basic test of which traces might be bad.
        [clustresults,CCCC,SUMDDD] = kmeans(all_datwf(sigwin,indgd)',2,'Distance','sqeuclidean','Replicates',20);
        num1 = length(find(clustresults == 1)); num2 = length(find(clustresults == 2));
        [minval,mindx] = min([num1 num2]); cvals = [1 2];
        mincluster= cvals(mindx);
        clustwf1 = CCCC(1,:);
        clustwf2 = CCCC(2,:);
        CMAT=corrcoef(clustwf1,clustwf2); 
        % relative phase is already obliterated in the mean stacks so shouldn't need to worry about correcting for lag
        % between the two cluster centroids
        if CMAT(1,2) > 0.65
            disp('Cluster Analysis is not useful right now.')
            donothing=1;
        else
            donothing=0;
        end
        
        % SNR calculation, windowing a period around the signal. 
        for isnr = 1:length(all_datwf(1,:))
           currwv = all_datwf(:,isnr);
           
           noisewin = find(-3*centperiod < ttws & ttws <  -1*centperiod); 
           sigwin   = find(-1*centperiod < ttws & ttws < 1.5*centperiod);
           
           SNR(isnr) = (rms(currwv(sigwin))/rms(currwv(noisewin)))^2;  
        end
        

        normf = ampf*max(max(abs(all_datwf(:,indgd))));
        % y option - distance or order
            
        if ydopt==1
            yds = [eqar(indgd).gcarc]; 
        else
            yds = 1:length(indgd);
        end
    
        yext = [min(yds),max(yds)];
        ax1 = axes(figxc,'pos',[0.05 0.07 0.3 0.85]); hold on
        cols = {'b','k'};
        for ig = 1:length(indgd)
            is = indgd(ig);
            currcluster=clustresults(ig);
            if currcluster == mincluster && donothing == 0 && clusteroff == 0
            % identify as the bad cluster. 
            plot(ax1,ttws,all_datwf(:,is)/normf/(10/diff(yext)) + yds(ig),'-','LineWidth',0.6,'color','r')
            else
            plot(ax1,ttws,all_datwf(:,is)/normf/(10/diff(yext)) + yds(ig),'-','LineWidth',0.6,'color',cols{isodd(ig)+1})
            end
            if exist('newxlims')
            text(ax1,newxlims(2)+1,yds(ig),datinfo(is).sta,'fontweight','bold','fontsize',20)                
            else
            text(ax1,prelimwind(2)+1,yds(ig),datinfo(is).sta,'fontweight','bold','fontsize',20)
            end
            if SNRdisp == 1 && exist('newxlims')
                 text(ax1,newxlims(1)+1,yds(ig),['SNR: ' num2str(SNR(is))],'fontsize',16)             
             
            elseif SNRdisp == 1
               text(ax1,prelimwind(1)+1,yds(ig),['SNR: ' num2str(SNR(is))],'fontsize',16)
             
            end

        end

        % if OverlayPhase 
        %     for phasechecker=1:length(PhasesforthisEQ)
        % 
        %         TmpATimelist = ArrivalTimeList(phasechecker,indgd);
        %         plot(TmpATimelist,yds,'color','k','linewidth',1)
        %     end
        % end

        title(ax1,sprintf('PICK XCOR WINDOW\nORID %u, seaz %.0f, gcarc %.0f\ndep %.0fkm, Ms %.1f',...
            orid,mean([eqar.seaz]),mean([eqar.gcarc]),evinfo.edeps(ie),evinfo.evmags(ie)),'Fontsize',18)
        if exist('newxlims')
        axis(ax1,[newxlims min(yds)-0.06*diff(yext) max(yds)+0.06*diff(yext) ])
        else
        axis(ax1,[prelimwind min(yds)-0.06*diff(yext) max(yds)+0.06*diff(yext) ])
        end
        set(ax1,'box','on','fontsize',16);
        try
            hp(1) = plot(ax1,[x1,x1],axlim(ax1,[3,4]),'--r','Linewidth',1.5);
            hp(2) = plot(ax1,[x2,x2],axlim(ax1,[3,4]),'--r','Linewidth',1.5);
        end

        %% Option to edit filter
        % work on figure 
        figure(figxc)
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
                if x2 < x1
                   disp('Wrong order for window start/endpoint selection')
                   try delete(hp); end
                   [x,junk]=ginput(1);
                    hp(1) = plot(ax1,[x1,x1],axlim(ax1,[3,4]),'--r','Linewidth',1.5);
                    [x2,y2] = ginput(1);
        
                end
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
                    
                    for jjkkllmm = 1:length(datinfo)
                        if isfield(datinfo,'xcor')
                        datinfo(jjkkllmm).xcor(FCounter) = 0;
                        else
                        datinfo(jjkkllmm).xcor = [0];                 
                        end
                    end
                    
                    %[datinfo.xcor] = deal(false);
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
                while isempty(flt)
                   disp('Enter new bounds for your filter! You have to press "ok"!')      
                                   flt = inputdlg({'Enter  max  period','Enter  min  period','Enter  N  filter  poles'},...
                               'Filter details',1,{num2str(1/filtfs(1)),num2str(1/filtfs(2)),num2str(npoles)});
                end
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
             case 't'  % display SNR

                if SNRdisp == 1
                 SNRdisp = 0;
                elseif SNRdisp == 0
                     SNRdisp = 1;
                
                end
                               
            case 'k' % toggle cluster analysis QC highlight switch on/off
                if clusteroff == 0
                    clusteroff = 1; %cluster analysis turned off now
                     disp('Cluster analysis toggled off!')
                   
                elseif clusteroff == 1
                    clusteroff = 0; %cluster analysis turned on now
                    disp('Cluster analysis toggled on!')
                end
                
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
            case 'c'
                input1 = input('Enter a new start time for the zoomed-in time window.');
                input2 = input('Enter a new end time for the zoomed-in time window.');
                while input2 < input1
                    disp('That''s the wrong order for the start/end times!')
                    input1 = input('Enter a new start time for the zoomed-in time window.');
                    input2 = input('Enter a new end time for the zoomed-in time window.');
                end
                newxlims=[input1 input2];
                
            case 'v' % make stacked spectogram  
                Make_Diagnostic_Figure
                
            case 'p' % overlay arrivals

                if OverlayPhase == 0
                         OverlayPhase=1;
                elseif  OverlayPhase==1       
                        OverlayPhase=0;
                end
             case 'z' % xlim anywhere
                k11=input('enter xlim start');
                 k12=input('enter xlim end');
               
                
        end


        % move on to xcorr if ENTER selected
        if isempty(b), cont = 1; end
    
    end % continue to process 
    
    % Salsipuede
	if b==27  % ESCape if escaped
        fprintf('Skipping event on user request...\n')
        for jjkkllmm = 1:length(datinfo)
        if isfield(datinfo,'xcor')
        datinfo(jjkkllmm).xcor = [datinfo(jjkkllmm).xcor 0];
        else
        datinfo(jjkkllmm).xcor = [0];                 
        end
        end
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
                        disp('did some deletions')
        elseif iter > 1
            indgd(acor<acormin) = []; % then kill low acor traces
            disp('did some deletions')
        end           
        if length(indgd) < 2,
            RESET = input('Move to next event? (0) Or Try again? (1)'); 
            break, 
        end
        iter=iter+1;
    end
	if length(indgd) < 2
        fprintf('NO GOOD TRACES/ARRIVALS, skip...\n')
        
        RESET = input('Move to next event? (0) Or Try again? (1)');
        if RESET
    load(arfile)      % loads eqar structure
   try cla(ax3); end
   try  cla(ax2); end
    try cla(ax1); end
    cont=0;
    % Best to once again decide which stations to delete for new sets of picks
    % at a different frequency. Sometimes bad ones aren't always bad. 
    indgd = 1:size(eqar); 
    indgd(mean(abs(all_dat0(:,indgd)))<1e-16)     = []; % kill zero traces
    indgd(mean(abs(all_datwf(:,indgd)))==0)     = []; % kill zero traces
    indgd(isnan(mean(abs(all_datwf(:,indgd))))) = []; % kill nan traces
    for ik = 1:length(manualkillstas),indgd(strcmp({eqar(indgd).sta},manualkillstas{ik})) = []; end % kill bad stas traces
    clear('x1','x2');
    clear newxlims

        else


                for jjkkllmm = 1:length(datinfo)
        if isfield(datinfo,'xcor')
        datinfo(jjkkllmm).xcor = [datinfo(jjkkllmm).xcor 0];
        else
        datinfo(jjkkllmm).xcor = [0];                 
        end
                end
                
        save(datinfofile,'datinfo');
        nextev = 1;
        continue

        end
    end
    if exist('RESET') == 0
        RESET = 0;
    end
    if RESET 
        %do nothing
    else
    %% ------------------ PLOTTING ------------------
    % RECORD SECTION TO CHECK
    try cla(ax3); end
    ax3 = axes(figxc,'pos',[0.68 0.07 0.27 0.85]); cla(ax3)
    set(ax3,'ytick',[],'box','on','fontsize',16,...
        'xlim',[ -25 45]); 
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
             case 'z' % Change the xlims here...
                k11=input('enter xlim start');
                 k12=input('enter xlim end');
               xlim([k11 k12])
        end
    end
	if isnan(conf)
                            for jjkkllmm = 1:length(datinfo)
                        if isfield(datinfo,'xcor')
                        datinfo(jjkkllmm).xcor(FCounter) = 0;
                        else
                        datinfo(jjkkllmm).xcor = [0];                 
                        end
                    end
        save(datinfofile,'datinfo');
        ie = ie + ied;
        nextev = 1;        
        continue
    end

    % while testing
    keyboard 
    return
    % while testing

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
    arfile2save      = [data_eqar_dir,evdir,'_EQAR_',phase,'_',component '_Filt' num2str(filtfs(1)) '_' num2str(filtfs(2)) '.mat'];
    
    % while testing
    keyboard 
    return
    % while testing

    save(arfile2save,'eqar')
    fprintf(' saved\n')
    % RECORD XCOR IN DATINFO 
    for jjkkllmm = 1:length(datinfo)
    if isfield(datinfo,'xcor')
    datinfo(jjkkllmm).xcor(FCounter) = 0;
    else
    datinfo(jjkkllmm).xcor= 0;                 
    end
    
    if isfield(datinfo,'xcor')
    datinfo(jjkkllmm).Freqlolist(FCounter) = filtfs(1);                 
    datinfo(jjkkllmm).Freqhilist(FCounter) = filtfs(2);       
    else
    datinfo(jjkkllmm).Freqlolist= filtfs(1);                 
    datinfo(jjkkllmm).Freqhilist= filtfs(2);                 
    end
    
    end
    
    for dxdx=indgd
    datinfo(dxdx).xcor(FCounter) = 1;
    end
    save(datinfofile,'datinfo')
    
    fprintf(' STA  CHAN  NEZ  resp  tilt  comp  xcor\n')
	for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f     %1.0f     %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).rmtilt,datinfo(is).rmcomp,datinfo(is).xcor); end
    
    
    
    %%% Check if we're done with picks or if we want to do them again at
    %%% another frequency. 
    Checker = input('Are you done with picks yet (Enter 1 at prompt), or would you like to make more measurements at another frequency (Enter 0 at prompt)?');
    if Checker
    % move on to next event
    ied = 1;
    ie = ie + ied;
    nextev = 1;
    ied = 1;
    try close figmetric; end
    else
            FCounter=FCounter+1;
    % load original EQAr and reset some configurations. 
    load(arfile)      % loads eqar structure
    figure(53)
    cla(ax3); 
     cla(ax2); 
     cla(ax1); 
    cont=0;
    % Best to once again decide which stations to delete for new sets of picks
    % at a different frequency. Sometimes bad ones aren't always bad. 
    indgd = 1:size(eqar); 
    indgd(mean(abs(all_dat0(:,indgd)))<1e-16)     = []; % kill zero traces
    indgd(mean(abs(all_datwf(:,indgd)))==0)     = []; % kill zero traces
    indgd(isnan(mean(abs(all_datwf(:,indgd))))) = []; % kill nan traces
    for ik = 1:length(manualkillstas),indgd(strcmp({eqar(indgd).sta},manualkillstas{ik})) = []; end % kill bad stas traces
    clear('x1','x2');
    clear newxlims
    end
     clf(figxc);

    end 

    try; clear RESET; end
    toc  
    end
end % while processing