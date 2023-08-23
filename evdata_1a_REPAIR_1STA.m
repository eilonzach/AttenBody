% script to go through data and reset to raw traces and channels
clear 

station = 'KIRE';
phases = {'P','S','SKS','PKS'};

overwrite = true;

ifsave = true;

% project details
dbname = 'EARdb';
dbdir = '/Users/zeilon/Work/BWAtten/EARdb/'; % include final slash

%% Preliminaries
wd = pwd;
addpath('matguts')
cd(dbdir);
run([dbdir,dbname,'_startup.m']);
javaaddpath('IRIS-WS-2.0.15.jar')
spd = 24*60*60; % seconds per day - to convert time in seconds to serial time

% event details
load([infodir,'/events'],'evinfo'); 
% station details
load([infodir,'/stations'],'stainfo');

orientation = 0;

for ie = 200:200 % 1:norids % loop on orids
    orid = evinfo.orids(ie);
    fprintf('\n Orid %.0f %s \n\n',orid,evinfo.evtimes_IRISstr{ie}(1:end-4))
    evdir = [num2str(orid,'%03d'),'_',evinfo.datestamp(ie,:),'/'];
    datinfofile = [datadir,evdir,'_datinfo'];

    load(datinfofile)
    if isempty(intersect({datinfo.sta}',station)), fprintf('This station has no data for this event\n'); continue; end % sta has no data for this evt
    is = find(strcmp({datinfo.sta}',station));
    
    load([datadir,evdir,station,'.mat']); 
    stadat = data;% load sta data for this evt
    chans = stadat.chans;
    raw = stadat.raw;
    
    %% do any corrections
    
    
    %% redo rotation
    ich = find(strcmp(chans.component,'H'));
    icn = find(strcmp(chans.component,'N'));if isempty(icn), icn = find(strcmp(chans.component,'1')); end
    ice = find(strcmp(chans.component,'E'));if isempty(ice), ice = find(strcmp(chans.component,'2')); end 
    icz = find(strcmp(chans.component,'Z'));

    ics = [ich icn ice icz];

    ichr = find(strcmp(raw.chans.component,'H'));
    icnr = find(strcmp(raw.chans.component,'N'));if isempty(icnr), icnr = find(strcmp(raw.chans.component,'1')); end
    icer = find(strcmp(raw.chans.component,'E'));if isempty(icer), icer = find(strcmp(raw.chans.component,'2')); end
    iczr = find(strcmp(raw.chans.component,'Z'));
    
    ic_new = [ichr icnr icer iczr];

    newdat = zeros(size(stadat.dat));
    sor = sind(orientation);
    cor = cosd(orientation);

    newdat(:,icn) = cor*raw.dat(:,icn) - sor*raw.dat(:,ice);    % NORTH
    newdat(:,ice) = sor*raw.dat(:,icn) + cor*raw.dat(:,ice);    % EAST
    newdat(:,ich) = raw.dat(:,ich);
    newdat(:,icz) = raw.dat(:,icz);

    newdat = newdat(:,ic_new);
    
    %% redo response removal
    for ic = 1:length(chans.component)
        fprintf(' %s,',chans.name{ic})
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



        %% RESPONSE REMOVAL, CRIBBED FROM JINGLE'S rm_resp SCRIPT
        lo_corner = 0.005;  % in Hz
        npoles=5;

        dat0 = newdat(:,ic);
        dat0 = detrend(dat0);
        dat0 = flat_hanning_win(1:length(dat0),dat0,1,length(dat0),50);

        T = stadat.tt(end)-stadat.tt(1);
        N = stadat.nsamps;
        delta = 1./stadat.samprate;
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

        lo_w=2*pi*lo_corner;
        hpfiltfrq=( ((w./lo_w).^(2*npoles))./(1+(w./lo_w).^(2*npoles)) );
        norm_trans=hpfiltfrq./resp;    % transfer function
        norm_trans(isnan(norm_trans)) = 0;

        fftdat0 = fft(dat0);
        fftdat0 = fftdat0(:).*norm_trans(:); % in f-space multiply (i.e. convolve time-domain function)
        dat0_cor = real(ifft(fftdat0));

        newdat(:,ic) = dat0_cor;
    end % loop on chans
        
        
    %% fix eqar structure
    stadat.dat    = newdat;
    stadat.NEZ    = true;
    stadat.rmresp = true;
    stadat.rmtilt = false;
    stadat.rmcomp = false;
    
    [stadat.phases.xtime]   = deal([]);
    [stadat.phases.xartime] = deal([]);
    [stadat.phases.xacor]   = deal([]);
    
    if ifsave
        data = stadat;
        save([datadir,evdir,station],'data')
    end
    
    %% fix datinfo
    datinfo(is).NEZ    = true;
    datinfo(is).rmresp = true;
    datinfo(is).rmtilt = false;
    datinfo(is).rmcomp = false;
    if ifsave
        save(datinfofile,'datinfo')
    end
    %% fix phase-specific structures
    fprintf('Phases: ')
    for ip = 1:length(phases)
        phase = phases{ip};
        fprintf('%s... ',phase)
        
        phdatinfofile = [datadir,evdir,'_datinfo_',phase,'.mat'];
        if isempty(dir(phdatinfofile)), continue, end
        load(phdatinfofile)
        is = find(strcmp({datinfo.sta}',station));
             
        eqarchans = {'Z','R','T'};
        for ic = 1:length(eqarchans)
            fprintf('%s',eqarchans{ic})
            eqarfile = [datadir,evdir,'_EQAR_',phase,'_',eqarchans{ic},'.mat'];
            load(eqarfile)

            if any(strcmp(stadat.chans.component,'Z'))
                eqar(is).datZ = interp1(stadat.tt,stadat.dat(:,strcmp(data.chans.component,'Z')),eqar(is).tt);
            end
            if any(strcmp(data.chans.component,'N')) & any(strcmp(data.chans.component,'E'))
                datN = interp1(stadat.tt,stadat.dat(:,strcmp(stadat.chans.component,'N')),eqar(is).tt);
                datE = interp1(stadat.tt,stadat.dat(:,strcmp(stadat.chans.component,'E')),eqar(is).tt);

                foraz = mod(stadat.seaz+180,360);
                eqar(is).datR =  datN*sind(foraz) + datE*sind(foraz);
                eqar(is).datT = -datN*sind(foraz) + datE*cosd(foraz);
            end
            
            if isfield(datinfo,'xcor')
            if datinfo(is).xcor
                eqar(is).dT = [];
                eqar(is).acor = [];
                eqar(is).abs_arrT = [];
                eqar(is).par_dT = [];
            end
            end
            if isfield(datinfo,'dtstar')
            if datinfo(is).dtstar
                eqar(is).snr_wf = [];
                eqar(is).frq = [];
                eqar(is).specn = [];
                eqar(is).specs = [];
                eqar(is).specss = [];
                eqar(is).fcross = [];
                eqar(is).dtstar = [];
                eqar(is).std_dtstar = [];
                eqar(is).par_dtstar_specR = [];
            end
            end
            
            if ifsave
                save(eqarfile,'eqar');
            end

        end % loop on eqarchans
        
        datinfo(is).NEZ    = true;
        datinfo(is).rmresp = true;
        datinfo(is).rmtilt = false;
        datinfo(is).rmcomp = false;
        datinfo(is).xcor = false;
        datinfo(is).dtstar = false;
        fprintf(' > ')
        if ifsave
            save(phdatinfofile,'datinfo');
        end
        
    end % loop on phases
    
    fprintf('fixed\n')
end % loop on events
cd(wd)