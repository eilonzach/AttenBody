% cycle through events and calculate each station's average relative
% amplification by looking at the spectral ratios, taking the average
% zero-frequency crossing points, and solving for each station's average
% zero-freq relative amplification value across all events
% 
% close all

% project details
dbname = 'EARdb';
dbdir = '/Users/zeilon/Dropbox/Work/EARdb/'; % include final slash


%% parameters
phase = 'S';
component = 'T'; %'Z', 'R', or 'T'

ifsave     = false;
ifplot     = true; 

%% conditions
mag_min = 6.25; % skip events if below this mag
acor_min = 0.75; % skip events if xcor max below this acor
snr_min = 3; % skip result if data below this snr


% datproc parms
resamprate = 4 ; % new, common sample rate
filtfs = 1./[40 1]; % [flo fhi] = 1./[Tmax Tmin] in sec
taperx = 0.2;
datwind = [-160 165]; % window of data in eqar structure
specwind = [-15 15]; % [-prex postx] So both values are relative to arrival (prex is thus +ive)
snrmin = 10;
spec2use = 'specss'; % 'specss' for smoothed spectrum, or 'specs' for unsmoothed
fmin = 0.045; % caps maximum value otherwise set by fcrosslo
fmax = 0.1;    % caps minimum value otherwise set by fcrosshi

firstorid = 400; % 412 is a good one for the amp2phiwt working, 2200


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

%% prep inversion structures
% model parameter will be list of stations, then list of events
M = stainfo.nstas + evinfo.norids;
si = [];
sj = [];
sv = [];
dlgA = [];
countGrows = 0;


for ie = 1:evinfo.norids  % loop on orids
    fprintf('Orid %.0f\n',ie)
    if  evinfo.evmags(ie)<mag_min, continue, end
    orid = evinfo.orids(ie);
    % name files and directories
    evdir = [num2str(orid,'%03d'),'_',evinfo.datestamp(ie,:),'/'];
    datinfofile = [datadir,evdir,'_datinfo_',phase];
    arfile      = [datadir,evdir,'_EQAR_',phase,'_',component];
       % check files exist
    if ~exist([datinfofile,'.mat'],'file'), fprintf('No station mat files for this event\n');continue, end
    if ~exist([arfile,'.mat'],'file'), fprintf('No arfile for this event + phase\n');continue, end    
    % load info file
    ld = load(datinfofile,'datinfo'); datinfo = ld.datinfo; % loads datinfo stucture
    
    % options to skip
    if isempty(datinfo), fprintf('No station mat files for this event\n'); continue, end
    if ~isfield(datinfo,'xcor'), fprintf('   NEED TO XCOR\n',phase), continue, end

    % load data file
    ld = load(arfile); eqar = ld.eqar;      % loads eqar structure
    if all([eqar.pred_arrT]==0), fprintf('No %s arrivals for this event\n',phase), continue, end
    
        
    %% parse stations that exist and have good spectra
    %  get mtm spec ratio data
    if ~isfield(eqar,'dtstar'), fprintf('   NEED TO CALC DTSTAR\n',phase), continue, end

    is_evt = []; % station number for this event - for indexing eqar
    is_db = []; % station number in entire database - for indexing model vector

    for is = 1:stainfo.nstas
        ise = find( strcmp({datinfo.sta},stainfo.stas(is)) & strcmp({datinfo.nwk},stainfo.nwk(is))  ); 
        if isempty(ise), continue; end
        if isempty(eqar(ise).specs), continue; end
        if eqar(ise).snr_wf < snr_min, continue; end
        if isempty(eqar(ise).dtstar), continue, end
        is_evt = [is_evt;ise];
        is_db = [is_db;is];
    end
    Nse = length(is_evt);

    % sta-by-sta acceptable f range from crossing freqs
    fcross = [[eqar(is_evt).fcrosslo]',[eqar(is_evt).fcrosshi]'];
    fcross(fcross(:,1) < fmin,1) = fmin;
    fcross(fcross(:,2) > fmax,2) = fmax;

%     figure(24),clf
%     subplot(211)
%     plot(eqar(is_evt(1)).frq,detrend(log([eqar(is_evt).specn]),0))
%     set(gca,'xscale','log')
%     subplot(212)
%     plot(eqar(is_evt(1)).frq,detrend(log([eqar(is_evt).(spec2use)]),0))
%     set(gca,'xscale','log')
%     pause


    % loop through pairs of stations, getting pairwise amplifications
    for is1 = 1:Nse
    for is2 = is1+1:Nse 
        spec1 = eqar(is_evt(is1)).(spec2use);
        spec2 = eqar(is_evt(is2)).(spec2use);
        lgspecR = log(spec2./spec1);
        
        % get frequency range for the spectral fit
        if ~isequal(eqar(is_evt(is1)).frq ,eqar(is_evt(is2)).frq) 
            warning('frequency vectors do not match between stations!'); 
        else
            frq = eqar(is_evt(is1)).frq;
        end % frequency vectors
        okfrq = (frq >= fcross(is1,1)) & (frq >= fcross(is2,1)) & (frq <= fcross(is1,2)) & (frq <= fcross(is2,2));
        lgspecR = lgspecR(okfrq);
        frq = frq(okfrq);
        Nf = sum(okfrq);
        
        % if not enough frequencies to estimate via straight line fit, just
        % take longest period spectral ratio
        if Nf<=2
            continue % or just skip!
            dAij = lgspecR(frq == min(frq)); 
        end 
        
        % straight line fit to get zero crossing
        Gsl = [ones(Nf,1),frq];
        msl = (Gsl'*Gsl)\Gsl'*lgspecR;
        dAij = msl(1);

        % some QC
        %  1) no negative dA! 
        %  2) no dA that's wildly different from the lowest freq dA
        % if dAij < 0 
        %     continue
        % else
            if abs(dAij - lgspecR(frq == min(frq))) > 1
            continue
        end

        %% insert into large sparse matrix index vectors
        countGrows = countGrows + 1;
        si = cat(1,si,countGrows*[1;1;1]);
        sj = cat(1,sj,[is_db(is1);is_db(is2); ie + stainfo.nstas]); % val = is2 - is1 + event term
        % actually don't want/need event term, as the measurements are
        % directly station to station
        sv = cat(1,sv,[-1;1;0]);
        dlgA(countGrows,1) = (dAij);

    end
    end
% if ie > 200, break; end
end % loop on orids

% build sparse G matrix
G = sparse(si,sj,sv,countGrows,M,3*countGrows);

nobspersta = full(sum(abs(G(:,1:stainfo.nstas)),1)');
nobsperevt = full(sum(abs(G(:,stainfo.nstas + [1:evinfo.norids]   )),1)');

% add constraint equation
G(countGrows+1,1:stainfo.nstas) = 1;
dlgA(countGrows+1,1)=0;


% solve!

% tiny bit of damping
mm = (G'*G + 1e-9*eye(M))\G'*dlgA;

% mm = lsqr(G,dlgA,1e-6,500);
dlgA_stas = mm(1:stainfo.nstas);
dlgA_evts = mm(stainfo.nstas+1:end);

% nan out ones with no information
dlgA_stas(nobspersta==0) = nan;
dlgA_evts(nobsperevt==0) = nan;


ev_i = sj(3:3:end) - stainfo.nstas;
s1_i = sj(1:3:end);
s2_i = sj(2:3:end);

% mean misfit by station
res = abs(dlgA(1:end-1)-G(1:end-1,:)*mm);
mean_res_sta = zeros(stainfo.nstas,1);
for is = 1:stainfo.nstas
    mean_res_sta(is) = mean(res(s1_i==is | s2_i==is));
end


%% -------------------------- PLOTS ---------------------------
if ifplot
    figure(9); clf, 
    scatter(stainfo.slons,stainfo.slats,10*sqrt(nobspersta)+0.01,dlgA_stas,'filled')
    caxis(log([0.2 5])), hc1 = colorbar;
    hcx = linspace(log(0.2),log(5),7)';
    set(hc1,'ytick',hcx,'yticklabel',cat(2,repmat('\times',7,1),num2str(round(exp(hcx),2))),'fontsize',15)
    title('station-wise average zero-f amplification (log scale)','fontsize',18)

    figure(10); clf, 
    scatter(evinfo.elons,evinfo.elats,10*sqrt(nobsperevt)+0.01,dlgA_evts,'filled')
    caxis(log([0.5 2])), hc2 = colorbar;
    hcx = linspace(log(0.5),log(2),7)';
    set(hc2,'ytick',hcx,'yticklabel',cat(2,repmat('\times',7,1),num2str(round(exp(hcx),2))),'fontsize',15)
    title('event-wise amplification correction (log_e units)','fontsize',18)

    figure(11);clf
    scatter(dlgA(1:end-1),G(1:end-1,:)*mm,20,res,'filled')
    xlabel('dlgA_{obs}','fontsize',16)
    ylabel('dlgA_{pre}','fontsize',16)
    axis equal
%     scatter(dlgA(1:end-1),abs(G(1:end-1,:)*mm - dlgA(1:end-1)),20,ev_i,'filled')

    figure(15);clf
    scatter(stainfo.slons,stainfo.slats,10*sqrt(nobspersta)+0.01,mean_res_sta,'filled')
    caxis([0 1]), colorbar
    title('station-wise average obs-pre |misfit|','fontsize',18)

end % ifplot

%% -------------------------- test output ---------------------------
unq_orid = unique(ev_i);
for ie = 1:length(unq_orid)  % loop on orids
    orid = evinfo.orids(unq_orid(ie));
    % name files and directories
    evdir = [num2str(orid,'%03d'),'_',evinfo.datestamp(unq_orid(ie),:),'/'];
    arfile      = [datadir,evdir,'_EQAR_',phase,'_',component];
    % load data file
    ld = load(arfile); eqar = ld.eqar;      % loads eqar structure
    Npaire = sum(ev_i == unq_orid(ie));
    if Npaire < handshake(4), continue; end % ignore evts with few stas

    figure(12);clf, hold on
    figure(13);clf, hold on
    ips = find(ev_i == unq_orid(ie));
    for ii = 1:length(ips) 
        
        ip = ips(ii);
        spec1 = eqar(strcmp({eqar.sta}',stainfo.stas(s1_i(ip)))).(spec2use);
        spec2 = eqar(strcmp({eqar.sta}',stainfo.stas(s2_i(ip)))).(spec2use);
        lgspecR = log(spec2./spec1);
        frq = eqar(strcmp({eqar.sta}',stainfo.stas(s1_i(ip)))).frq;
        A0_obs = dlgA(ip);
        A0_est = G(ip,:)*mm;
        figure(12)
        plot(frq,lgspecR,'-xk')
        plot(0.0001,A0_est,'ob')
        plot(0.0001,A0_obs,'or')
        xlim([0 0.2])
        figure(13)
        plot(A0_obs,A0_est,'ok')
        text(exp(dlgA(ip))+0.01,exp(G(ip,:)*mm)+0.01,sprintf('%.0f,%.0f',s1_i(ip),s2_i(ip)))
        pause
    end
end




%% -------------------------- SAVE ---------------------------
if ifsave
    amplif = struct('stas',{stainfo.stas},'nwk',{stainfo.nwk},...
                'dlgA_stas',dlgA_stas,'nobs_dlgA_s',nobspersta,'meanL1misfit',mean_res_sta);
%     amplif = struct('stas',{stainfo.stas},'nwk',{stainfo.nwk},...
%                 'dlgA_stas',dlgA_stas,'nobs_dlgA_s',nobspersta,...
%                     'orid',evinfo.orids,...
%                 'dlgA_evts',dlgA_evts,'nobs_dlgA_e',nobsperevt);
    ofile = sprintf('%sA0amplification_%s%s',resdir,phase,component);
    save(ofile,'amplif')

end 

cd '/Users/zeilon/Dropbox/MATLAB/AttenBody';
