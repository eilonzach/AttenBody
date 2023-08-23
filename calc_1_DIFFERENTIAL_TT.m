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
%       y	toggle between gcarc distance and simply station order
%       d	delete station with trace beneath cursor
%       u	flip sign of station with trace beneath cursor
%    left   prev event
%   right 	next event 
%      up 	increase amplitudes
%    down 	decrease amplitudes
%     ESC   stop

clear

% % project details
% dbname = 'EARdb';
% dbdir = '/Users/zeilon/Work/EastAfrica/EARdb/'; % include final slash

% project details
dbname = 'FRES_PILOT';
dbdir = '~/Dropbox/Work/FRES_PILOT/'; % include final slash


%% parameters
phase = 'SKS';
component = 'R'; %'Z', 'R', or 'T'
resamprate = 10 ; % new, common sample rate
% Do filtfs [12 1] for S, and [5 1] for P
% filtfs = 1./[2 0.25]; % P [flo fhi] = 1./[Tmax Tmin] in sec 
filtfs = 1./[12 1]; % S [flo fhi] = 1./[Tmax Tmin] in sec 
taperx = 0.2;
npoles = 2;
datwind = [-200 200]; % window of data in eqar structure
prelimwind = [-40 40];
% prelimwind = [-25 15];
xcorlagmax = 20;
gcdistmax = 12;
mingcarc = 30;
maxgcarc = 130;
acormin = 0.75;
overwrite = false;
minNtrace = 4;
% FRESPILOT -  S(T) got to 1863. P(Z) got to 1920. SKS(R) got to 1847
% EARdb (I think outdated) SKS got to 1044; S got to 1287; % need to do 160-335 again on all
firstev = 1847; 


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

ie = firstev; % up to 364
ied = 1; %incrementer - leave as 1
while 1 %evinfo.norids % loop on orids 
%     if  mags(ie)<6.9, continue, end
    tic
    close all, clear('x1','x2');
    orid = evinfo.orids(ie);
    fprintf('\n Orid %.0f %s \n\n',orid,evinfo.evtimes_IRISstr{ie}(1:end-4))
    % name files and directories
    evdir = [num2str(orid,'%03d'),'_',evinfo.datestamp(ie,:),'/'];
    datinfofile = [data_eqar_dir,evdir,'_datinfo_',phase];
    arfile      = [data_eqar_dir,evdir,'_EQAR_',phase,'_',component];
   
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
    indgd(mean(abs(all_dat0(:,indgd)))<1e-16)     = []; % kill zero traces
    indgd(mean(abs(all_datwf(:,indgd)))==0)     = []; % kill zero traces
    indgd(isnan(mean(abs(all_datwf(:,indgd))))) = []; % kill nan traces
    for ik = 1:length(manualkillstas),indgd(strcmp({eqar(indgd).sta},manualkillstas{ik})) = []; end % kill bad stas traces
    
    if strcmp(phase,'SKS') ,indgd([eqar(indgd).gcarc]<86) = []; end % kill close traces
    
    if length(indgd) < minNtrace, fprintf('NOT ENOUGH GOOD TRACES/ARRIVALS, skip...\n'); ie = ie + ied; continue, end
    
    fprintf('Orid %.0f %s \n',orid,evinfo.evtimes_IRISstr{ie}(1:end-4))
	%% tauptime to help picking
    tauptime('d',mean([eqar.gcarc]),'z',evinfo.edeps(ie))
  
%% ------------------ PICK XCOR WINDOW, RECLEAN ------------------
    figxc = figure(53); clf(figxc); set(figxc,'position',[710 1 1626 1344])
    
    %% PICK TRAVEL TIMES
    cont=0;
    nextev=0;
    ydopt = 1;
    ampf = 1;
    while nextev==0
    while cont==0
    
        try delete(ax1); end


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
            plot(ax1,ttws,all_datwf(:,is)/normf/(10/diff(yext)) + yds(ig),'-','LineWidth',0.6,'color',cols{isodd(ig)+1})
            text(ax1,prelimwind(2)+1,yds(ig),datinfo(is).sta,'fontweight','bold','fontsize',20)
        end
        title(ax1,sprintf('PICK XCOR WINDOW\nORID %u, seaz %.0f, gcarc %.0f\ndep %.0fkm, Ms %.1f',...
            orid,mean([eqar.seaz]),mean([eqar.gcarc]),evinfo.edeps(ie),evinfo.evmags(ie)),'Fontsize',18)
        axis(ax1,[prelimwind min(yds)-0.06*diff(yext) max(yds)+0.06*diff(yext) ])
        set(ax1,'box','on','fontsize',16);
        try
            hp(1) = plot(ax1,[x1,x1],axlim(ax1,[3,4]),'--r','Linewidth',1.5);
            hp(2) = plot(ax1,[x2,x2],axlim(ax1,[3,4]),'--r','Linewidth',1.5);
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
        'xlim',[ -10 30]); 
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