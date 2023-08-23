% script to go through data and correct for instrument response for  each 
% of the station's channels. 
% 
% ONE APPROACH - SACPZ is to get instrument response files in SAC poles and
% zeros format, with file naming style:
%   SAC_PZs_NWK_STA_CHAN 
% where NWK is the 2-character network name, STA is the name of the
% station, and CHAN is the channel code.
% 
% if using rdseed, try:
% rdseed -p -f DatalessSEEDfilename
%
% AUTOMATICALLY CORRECTS TO VELOCITY (removing a zero from the IRIS SAC_PZ)
% 
% OTHER APPROACH - fetchPZ is to get instrument response along with the
% irisFetch station/channel request.
% 
% Z. Eilon 2016, updated 03/2023
 
clear 

resp_opt = 'fetchPZ'; % 'fetchPZ' or 'SACPZ'
overwrite = true;
ifplot = 0; 
startorid = 1151;

% % project details
% dbname = 'EARdb';
% dbdir = '/Volumes/Lacie/Granite_EastAfrica/EARdb/'; % include final slash

% project details
dbname = 'FRES_PILOT';
dbdir = '~/Dropbox/Work/FRES_PILOT/'; % include final slash


%% Preliminaries
wd = pwd;
addpath('matguts')
cd(dbdir);
run([dbdir,dbname,'_startup.m']);
spd = 24*60*60; % seconds per day - to convert time in seconds to serial time

% point to right data dir 
if exist(datadir,'dir')~=7, datadir = regexprep(datadir,'data','data-1'); end
% event details
load([infodir,'/events'],'evinfo'); 
% station details
load([infodir,'/stations'],'stainfo');

for ie = startorid:evinfo.norids %  loop on orids
    orid = evinfo.orids(ie);
    fprintf('\n Orid %.0f %s \n\n',orid,evinfo.evtimes_IRISstr{ie}(1:end-4))
    evdir = [num2str(orid,'%03d'),'_',char(evinfo.datestamp(ie,:)),'/'];
    datinfofile = [datadir,evdir,'_datinfo'];
   
    % check files exist
    if ~exist([datinfofile,'.mat'],'file'), fprintf('No data at all for this event\n');continue, end
    load(datinfofile)
    if isempty(datinfo), fprintf('No station mat files for this event\n');continue, end

    for is = 1:length(datinfo) % loop on stas
        sta = datinfo(is).sta; % sta name
        fprintf('Station %.0f %-5s...',is,sta)
        if datinfo(is).rmresp && ~overwrite % skip if already removed response
            fprintf(' done already\n'), continue, 
        end 

%         if ~iscell(respfiles), respfiles = {respfiles}; end
%         if ~iscell(respdirs), respdirs = {respdirs}; end
        
        load([datadir,evdir,sta,'.mat']); % load sta data for this evt
        if data.rmresp && ~overwrite % skip if already removed response
            fprintf(' done already\n'), continue, 
        end 
        
        chans = data.chans;
        chans_raw = data.raw.chans;
        
        ich = find(strcmp(chans.component,'H'));
        icn = find(strcmp(chans.component,'N'));
        ice = find(strcmp(chans.component,'E')); 
        icz = find(strcmp(chans.component,'Z'));
        if isempty(icn), icn = find(strcmp(chans.component,'1')); end
        if isempty(ice), ice = find(strcmp(chans.component,'2')); end
        ics = [ich icn ice icz];
        
        ichr = find(strcmp(chans_raw.component,'H'));
        icnr = find(strcmp(chans_raw.component,'N'));
        icer = find(strcmp(chans_raw.component,'E')); 
        iczr = find(strcmp(chans_raw.component,'Z'));
        if isempty(icnr), icnr = find(strcmp(chans_raw.component,'1')); end
        if isempty(icer), icer = find(strcmp(chans_raw.component,'2')); end
        icsr = [ichr icnr icer iczr];
        
        for ic = 1:length(chans.component)
            fprintf(' %s,',chans.name{ic})
            % icsr == ics(ic) % Need this chicanery in case chnames changed


            switch resp_opt 
                case 'SACPZ'
                    % find resp file
                    respfile_pre = ['SAC_PZs_',datinfo(is).nwk,'_',datinfo(is).sta,'_',chans.name{ic}];
                    respfile = dir([respdir,respfile_pre,'*']);
                    if length(respfile) > 1 % multiple response files for this nwk+sta+chan
                        for ir = 1:length(respfile)
                        % get beginning and end of response file in serial time
                        [~,~,~,~,~,resp_begin,resp_end]=parse_sacpz_filename(respfile(ir).name);
                        [rbm,rbd] = calday(resp_begin(1),resp_begin(2));
                        [rem,red] = calday(resp_end(1),resp_end(2));
                        rb = datenum([resp_begin(1),rbm,rbd,resp_begin(3:end)]);
                        re = datenum([resp_end(1),rem,red,resp_end(3:end)]);
                        if (rb < evinfo.evtimes(ie)) && (re > evinfo.evtimes(ie))
                            % if event within this time, this is desired respfile
                            break;
                        end
                        end
                    elseif length(respfile)==1
                        ir = 1;
                    else
                        fprintf('NO RESP'), continue
                    end
                    
                    % Helen's function
                    % [zz,pp,gain] = read_sac_pole_zero([respdir,respfile(ir).name]);
                   
                    % use function from seizmo toolbox
                    [pz]=readsacpz_rdseed([respdir,respfile(ir).name]);
                    gain = pz.k;
                    zz = complex(pz.z{1});
                    pp = complex(pz.p{1});
                    ounit = pz.input{1};
                    if strcmp(ounit,'M') % PZ file set up to convert to displacement
                        % correct to VELOCITY by removing a zero
                        zz = zz(2:end);
                    end     
                case 'fetchPZ'
                    %% pull out relevant info from stainfo
                    % time is evinfo.evtimes(ie)
                    % chan is chans.name{ic}
                    % station is sta
                    % network is datinfo(is).nwk

                    % test if the info is there in the stainfo structure
                    if ~isfield(stainfo,'chanresps')
                        error('No responses in this stainfo file... cannot use fetchPZ approach')
                    end
                    % find right station in stainfo
                    ista = find( strcmp(stainfo.stas,datinfo(is).sta) & strcmp(stainfo.nwk,datinfo(is).nwk) );
                    if length(ista) > 1
                        keyboard
                    end
                    ichn = find(strcmp(stainfo.chans(ista,:),chans.name{ic}));
                    % this is the response for this station's channel
                    responses = stainfo.chanresps{ista,ichn};
                    % now find correct time
                    ires = (evinfo.evtimes(ie) > [responses.StartDate_ser]) & ...
                           (evinfo.evtimes(ie) < [responses.EndDate_ser]);
                    if ~any(ires), fprintf('NO RESP'), continue; end
                    
                    response = responses(ires);
                    gain = response.Gain;
                    zz = response.Zeros;
                    pp = response.Poles;
                    ounit = response.outUnits;
                    if strcmp(ounit,'M') % PZ file set up to convert to displacement
                        % correct to VELOCITY by removing a zero
                        zz = zz(2:end);
                    end     
                    % while testing
%                     plot_resp_pz(pp,zz,gain)
                end

            
            %% RESPONSE REMOVAL, CRIBBED FROM JINGLE'S rm_resp SCRIPT
            % options
            lo_corner = 0.005;  % in Hz
            npoles=5;
            
            dat0 = data.dat(:,ic);
            dat0 = detrend(dat0);
            dat0 = flat_hanning_win(1:length(dat0),dat0,1,length(dat0),50);
            
            isserialtime = data.tt(1)>7e5;
            T = data.tt(end)-data.tt(1);
            if isserialtime, T = T*24*3600; end % need seconds!!!
            N = data.nsamps;
            delta = 1./data.samprate;
            if mod(N,2)
                 faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/T);
            else
                 faxis = [0:N/2,-N/2+1:-1]*(1/T);
            end
            w = faxis.*2*pi;
            resp = ones(size(w));
            for ip = 1:length(pp)
                resp = resp./(1i*w - pp(ip));
            end
            for ip = 1:length(zz)
                resp = resp.*(1i*w - zz(ip));
            end
            resp = resp*abs(gain); % made abs to get rid of neg gain issues!

            if ifplot
                figure(33)
                clf
                set(gcf,'position',[360   514   900   400]);
                hold on
                subplot(1,2,1)
                set(gca,'fontsize',18)
                semilogy(faxis,abs(resp),'rx');
                subplot(1,2,2)
                set(gca,'fontsize',18)
                plot(faxis,angle(resp),'rx');
            end
            
            lo_w=2*pi*lo_corner;
            hpfiltfrq=( ((w./lo_w).^(2*npoles))./(1+(w./lo_w).^(2*npoles)) );
            norm_trans=hpfiltfrq./resp;    % transfer function
            norm_trans(isnan(norm_trans)) = 0;

            fftdat0 = fft(dat0);
            fftdat0 = fftdat0(:).*norm_trans(:); % in f-space multiply (i.e. convolve time-domain function)
            dat0_cor = real(ifft(fftdat0));
            
            if ifplot
                figure(2); clf; hold on
%                 plot( data.dat(:,ic)./max(abs(data.dat(:,ic))) - 1,'b')
%                 plot( dat0./max(abs(dat0)),'r')
%                 plot( dat0_cor./max(abs(dat0_cor)) + 1,'g')
                d1 = data.dat(:,ic); %d1 = filt_quick(d1,0.01,.2,1./data.samprate);
                d2 = dat0;           %d2 = filt_quick(d2,0.01,.2,1./data.samprate);
                d3 = dat0_cor;       %d3 = filt_quick(d3,0.01,.2,1./data.samprate);
                d4 = data.raw.dat(:,ic);%d3 = filt_quick(d3,0.01,.2,1./data.samprate);
                plot( d1./max(abs(d1)) - 1,'b')
                plot( d4./max(abs(d4)) - 2,'k')
                plot( d2/max(abs(d2)),'r')
                plot( d3/max(abs(d3)) + 1,'g')
            end
            
            % while debugging 
            if ifplot
                figure(144); clf
                subplot(2,1,1);plot(data.dat(:,ic))
                subplot(2,1,2);plot(dat0_cor)
            end


            data.dat(:,ic) = dat0_cor;
        end % loop on chans
        
        
        %% log resp removal and save
        fprintf(' resp removed\n')
        data.rmresp = true;
        save([datadir,evdir,sta],'data')
        
        datinfo(is).rmresp = true;
        save(datinfofile,'datinfo')
                    
    end % loop on stas

	%% sum up
    fprintf(' STA  CHAN  NEZ  resp  tilt  comp\n')
    for is = 1:length(datinfo), fprintf('%-5s %4s   %1.0f     %1.0f     %1.0f     %1.0f\n',datinfo(is).sta,[datinfo(is).chans{:}],datinfo(is).NEZ,datinfo(is).rmresp,datinfo(is).rmtilt,datinfo(is).rmcomp); end
            
end
cd(wd);
