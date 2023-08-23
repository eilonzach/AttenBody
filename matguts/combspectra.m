function [delta_tstar,delta_T,std_dtstar,pairwise,fmids,delta_A0] = combspectra(dat,fcross,samprate,parms,ifplot,mtmspecss)
% [delta_tstar,delta_T,std_dtstar,pairwise,fmids,delta_A0] = combspectra(dat,fcross,samprate,parms,ifplot,mtmspecss)

%% establish some parms
if nargin < 5 
    ifplot = false;
end
if nargin < 6
    parms.qc.mtmRcomp = 0;
    mtmspecss = [];
end

Tmin = parms.comb.Tmin;
Tmax = parms.comb.Tmax;
Nwds = parms.comb.Nwds;
Tw_opt = parms.comb.Tw_opt;
npol = parms.comb.npol;

pretime = parms.wind.pretime;
prex    = parms.wind.prex;
postx   = parms.wind.postx;
taperx  = parms.wind.taperx;

minacor = parms.qc.minacor;
maxphi  = parms.qc.maxphi;
try
    maxdA = parms.qc.maxdA;
catch
    maxdA = 0.3;
end

amp2phiwt = parms.inv.amp2phiwt;
fmin_cskip= parms.inv.fmin_cskip;
% fmax      = parms.inv.fmax;
corc_skip = parms.inv.corr_c_skip;
ifwt      = parms.inv.ifwt;
if isfield(parms.inv,'R2default')
    R2default = parms.inv.R2default;
else
    R2default = 0; % if not set, ignore. 
end

if ~isfield(parms.inv,'alpha')
    alpha = 0;
elseif isempty(parms.inv.alpha)
    alpha = 0;
else
    alpha = parms.inv.alpha;
end

if ~isfield(parms.qc,'mtmRcomp') % did we specify whether to compare amplitudes to specR
    parms.qc.mtmRcomp = false;
elseif isempty(parms.qc.mtmRcomp)
    parms.qc.mtmRcomp = false;
end

fnq = samprate/2;
dt = 1./samprate;
npt = size(dat,1);

Nstas = size(dat,2);
if Nstas == 1, error('Only 1 station''s data!'), end
N = handshake(Nstas);

dtstar_pairwise = zeros(N,1);
dT_pairwise = zeros(N,1);
dlgA0_pairwise = zeros(N,1);
misfitnormed_pairwise = zeros(N,1);
R2_pairwise = zeros(N,1);
dtstar_std_pairwise = zeros(N,1);
dT_std_pairwise = zeros(N,1);
dlgA0_std_pairwise = zeros(N,1);

As_pairwise = zeros(N,Nwds);
phis_pairwise = zeros(N,Nwds);
wts_pairwise = zeros(N,Nwds);
inds_pairwise = zeros(N,Nwds);

%% prepare filter + cleaning parms
% Make set of period windows for bandpass filter
Tmids = logspace(log10(Tmin),log10(Tmax),Nwds)';
if strcmp(Tw_opt,'scale')
%%  ZE EDIT 5/18/22 - WAS
% Twdhs = 0.5*diff(logspace(log10(Tmin/2),log10(2*Tmax),Nwds+1)');
    %%  ZE EDIT 5/18/22 - NOW
    Twdhs = 0.5*diff(logspace(log10(Tmin/2),log10(Tmax/2),Nwds+1)');
else
    Twdhs = 0.5*Tw_opt(:).*ones(size(Tmids));
end


% Twdhs = 0.5*Tmids;
fmids = 1./Tmids;

flos = 1./(Tmids + Twdhs);
fhis = 1./(Tmids - Twdhs);

%% plot filters
if ifplot~=0
    plot_filter_comb( flos,fhis,dt,npol )
end

% calculate bandpass filter parameters
clear fltinfo
for iw = 1:Nwds
    % option 1: butter
%     [fltdat(iw).bb,fltdat(iw).aa]=butter(parms.comb.npol, [flos(iw), fhis(iw)].*dt.*2);
    % option 2: cheby
%     [fltdat(iw).bb,fltdat(iw).aa]=cheby1(parms.comb.npol,0.5, [flos(iw), fhis(iw)].*dt.*2);
    % option 3: higher order butter
    [z,p,k]=butter(npol, [flos(iw), fhis(iw)].*dt.*2.);
    [sos,g]=zp2sos(z,p,k); 
    fltinfo(iw).bb=sos; 
    fltinfo(iw).aa=g;
    fltinfo(iw).fmid = fmids(iw);
    fltinfo(iw).flo = flos(iw);
    fltinfo(iw).fhi = fhis(iw);
end

% set up cleaning/taper/window parms
nwin=round((postx+prex)/dt); % window length in samples

wdo1 = tukeywin(npt,2*taperx); % tapered for the filtering of the whole signal

n1=round((pretime-prex)/dt); % first sample in window
wdo2=[zeros(n1,1);tukeywin(nwin,2*taperx); zeros(npt-n1-nwin,1)]; % taperx% tukey window
ibds=[max(1,n1-floor(nwin*4*taperx)), min(npt,n1+nwin+floor(nwin*4*taperx))]; % pad/extend so traces are extra 4xtaperF (often 40%) of window on either side
jbds=ibds(1):ibds(2); % indices of points to keep
    
% clean the data (shouldn't need it!)
nnan=~isnan(dat);
dat = reshape(detrend(dat(nnan)),size(dat)); % detrend non-nan data
dat(isnan(dat)) = 0; %set any nans to zero


%% loop over station pairs
ui = zeros(2*N,1);
uj = zeros(2*N,1);
u  = zeros(2*N,1);
count = 0;

% hw = waitbar(count/N,'Progress through station-station dA,dphi');

for is1 = 1:Nstas
for is2 = is1+1:Nstas
%     waitbar(count/N,hw)
    count = count+1;
    % make elements of eventual G matrix
    ui(2*(count-1)+[1 2]) = count;
    uj(2*(count-1)+[1 2]) = [is1 is2];
    u (2*(count-1)+[1 2]) = [-1 1]; % delta is value of 2 - value of 1
    
    %% Do the work of running through the comb - option to correct cycle skip
    [ As1,phis1,wts1 ] = run_comb( dat(:,is1),dat(:,is2),fltinfo,wdo1,wdo2,jbds,dt,pretime,maxphi,corc_skip,fmin_cskip,ifplot==2 );
    [ As2,phis2,wts2 ] = run_comb( dat(:,is2),dat(:,is1),fltinfo,wdo1,wdo2,jbds,dt,pretime,maxphi,corc_skip,fmin_cskip,0 );
    % convert As2 into its inverse (account for the difference in mean
    % values due to frequency-independent differences in gains...)
% %     As2_to_1 = 1./As2/mean([mean(As1),1./mean(As2)]);
%     As2_to_1 = mean(As1.*As2)./As2; 
    As2_to_1 = 1./As2;
    phis2_to_1 = -phis2; % phis2 is opposite in sign to phis1, so this is the mean 
    % combine for As(1 vs 2; i.e. like As1) with simple means
    As = (As1+As2_to_1)/2;
    dAs = abs((As1-As2_to_1)./As); % fractional difference in amplitude (should be zero in ideal case)
    phis = (phis1+phis2_to_1)/2; 
    wts = (wts1+wts2)/2;
    
    % plot difference in two directions
    if ifplot
        figure(1001),clf; 
        subplot(311);hold on;
        scatter(log(As1),log(As2_to_1),wts*60,'b','filled'); xlabel('ln(As_1)'),ylabel('ln(As_2)')
        scatter(log(As1),log(As),wts*60,'c','filled'); 
        plot([0,min([2,max([As1;As2_to_1])])],[0,min([2,max([As1;As2_to_1])])],'--r'); 
%         xlim([0,min([2,max([As1;As2_to_1])])]); ylim([0,min([2,max([As1;As2_to_1])])]);
        subplot(312);
        scatter(phis1,phis2_to_1,wts*60,'b','filled'); xlabel('phis_1'),ylabel('phis_2')
        subplot(313);
        scatter(wts1,wts2,60,'b','filled'); xlabel('wts_1'),ylabel('wts_2')
    end
    
    if ~ifwt, wts(iw) = 1; end
    wts(As<0) = 0;
    As(As<0) = 0;
    wts(phis==maxphi) = 0;
    
    % work out frequency limit for this pair
    fmin_pair = max(fcross([is1,is2],1));
    fmax_pair = min(fcross([is1,is2],2));

    % if asked, compare the measured amplitude values to the spectral ratio
    % equivalent
    if parms.qc.mtmRcomp
        As_mtmR = interp1(mtmspecss.frq,mtmspecss.specs(:,is2)./mtmspecss.specs(:,is1),fmids);
        wts_mtmRcomp = .5./(0.25 + abs(log(As) - log(As_mtmR))); % water level it to 0.5. 
        wts_mtmRcomp = wts_mtmRcomp(:);
        % if there is ZERO difference in log amplitudes, then gets double
        % weighted. If there is 0.25 ln unit difference, gets normal weight.
        % If there is 1.25 difference in weight, down to third weight
    end

    %% QC
    inds = (fmids>=fmin_pair) & ...
           (fmids<=fmax_pair) & ...
           (abs(phis)<maxphi) & ...
           (sqrt(wts)>minacor) & ...
           (dAs<maxdA);
       
    % add dA assessment to weights
    wts = wts.*(1 - dAs./(1.5*maxdA));
    % add diff from specR to weights
    if parms.qc.mtmRcomp
        wts = wts.*wts_mtmRcomp;
    end

    pwts = wts; % save for plotting
    % kill weights of pts that don't pass QC
	wts(~inds) = 0;

    
    if sum(inds)<4 % only do if there are at least 4 ok measurements!
        As_pairwise(count,:) = As;
        phis_pairwise(count,:) = phis;
        wts_pairwise(count,:) = zeros(size(wts));
        
        dtstar_pairwise(count) = nan;
        dT_pairwise(count) = nan;
        misfitnormed_pairwise(count) = Inf;
        continue % skip to next station pair
    end
    
	% re-norm wts
    wts(inds) = wts(inds)./rms(wts(inds)); % DANGER - if all wts zero, gives inf. Should not have all zero - see above if statement.
    if any(isinf(wts))
        keyboard
    end

    
    %% calculate dtstar and dT
    % estimate from simultaneous inversion of amp and phase data
    [ dtstar,dT,A0,misfit,resnorm,A_SSE,P_SSE,lnA_SSE,dtstar_std,dT_std,dlnA0_std ] ...
        = invert_1pair_Aphi_4_dtdtstar( As(inds),phis(inds),fmids(inds), wts(inds),amp2phiwt,alpha);
    
    % Calculate scaled goodness of fit (R2) to affect weights for inv
    % RMS = 1 - SSE/SSTotal (both weighted)
    R2_phi = 1 - P_SSE./(phis(inds)'*diag(wts(inds))*phis(inds));
    R2_lnA = 1 - lnA_SSE./(log(As(inds))'*diag(wts(inds))*log(As(inds)));
    R2 = (amp2phiwt*R2_lnA + R2_phi)/(amp2phiwt+1);
    if isinf(R2) || isnan(R2) || R2 < 0 % R2<0 means made fit worse than null fit of zeros
        R2=0;
    end
    if R2default == 0 
        R2fac = 1;
    elseif R2default > 0 && R2default <= 1
        R2fac = R2/R2default;
    else
        warning('inappropriate R2 default - should be 0 to ignore or 0<R2<1');
    end
    
    %% plots
    if ifplot~=0
        figure(4), clf, set(gcf,'pos',[600 600 500,700])

        subplot(211), hold on
        pwts(pwts<=0)=nan;
        scatter(fmids,log(As),110*pwts,'or','MarkerFaceColor','r')
        scatter(fmids(inds),log(As(inds)),110*pwts(inds),'ok','MarkerFaceColor','r')
%         scatter(fmids(dodo),log(As(dodo)),100*wts(dodo),'o','MarkerFaceColor','b')
        if alpha==0
            plot(fmids,log(A0) - pi*fmids*dtstar,'g','Linewidth',1.5)
        else
            plot(fmids,log(A0) - pi*fmids.^(1-alpha)*dtstar,'g','Linewidth',1.5)
        end
        plot(fmin_pair*[1 1],[-1 1],'--b')
        plot(fmax_pair*[1 1],[-1 1],'--b')
        xlabel('freq','FontSize',18), ylabel('log(Amp)','FontSize',18)
        title(sprintf('Station %.0f vs. station %.0f',is1,is2),'FontSize',20)
        set(gca,'Xscale','linear','xlim',[0.03 1])
        if ~isempty(mtmspecss)
            plot(mtmspecss.frq,log(mtmspecss.specss(:,is2)./mtmspecss.specss(:,is1)),'-xk')
        end

        subplot(212), hold on
        scatter(fmids,phis,110*pwts,'or','MarkerFaceColor','r')
        scatter(fmids(inds),phis(inds),110*pwts(inds),'ok','MarkerFaceColor','r')
%         scatter(fmids(dodo),phis(dodo),100*wts(dodo),'o','MarkerFaceColor','b')
        if alpha==0
            plot(fmids,(log(fnq) - log(fmids))*dtstar./pi + dT,'g','Linewidth',1.5)
        else
            plot(fmids,0.5 * cot(alpha*pi/2) * fmids.^-alpha * dtstar + dT,'g','Linewidth',1.5)
        end
        plot(fmin_pair*[1 1],[-1 1],'--b')
        plot(fmax_pair*[1 1],[-1 1],'--b')
        set(gca,'Xscale','log','xlim',[0.03 1.],'ylim',maxphi*[-1 1]/2)
        title(sprintf('Misfit/N = %.3f  tstar=%.3f  dT=%.3f',misfit./length(inds),dtstar,dT),'FontSize',20)
        xlabel('log(freq)','FontSize',18), ylabel('phase-shift (s)','FontSize',18)
    end
    
    As_pairwise(count,:) = As;
    phis_pairwise(count,:) = phis;
    wts_pairwise(count,:) = wts * R2fac;
    inds_pairwise(count,:) = inds;

    dtstar_pairwise(count) = dtstar;
    dT_pairwise(count) = dT;
    dlgA0_pairwise(count) = log(A0);
    misfitnormed_pairwise(count) = misfit./sum(inds)./R2fac; % weight  will be  1./misfit, normalised by number of datapoints, with optional R2 scaling
    R2_pairwise(count) = R2; % weighted, amp2phweighted R2 (1=perfect fit)

    dtstar_std_pairwise(count) = dtstar_std;
    dT_std_pairwise(count) = dT_std;
    dlgA0_std_pairwise(count) = dlnA0_std;


end
end
% delete(hw)

% make structure of all pair-wise comparisons (Amp and phi spectra, as well
% as dtstars,dTs from individual pairs etc.)
pairwise = struct('dtstar',dtstar_pairwise,'dT',dT_pairwise,'misfit_normed',misfitnormed_pairwise,'R2_pairwise',R2_pairwise,...
                  'As',As_pairwise,'phis',phis_pairwise,'wts',wts_pairwise,'inds',inds_pairwise);

%% solve the least squares problem
G = sparse(ui,uj,u,N,Nstas,2*N);
    
covd_ts = diag(dtstar_std_pairwise.^2);
covd_T = diag(dT_std_pairwise.^2);
covd_A = diag(dlgA0_std_pairwise.^2);

if ifwt
%     W = diag(1./misfitnormed_pairwise); % old way
    W_ts = diag(diag(covd_ts).^-1);
    W_T  = diag(diag(covd_T).^-1);
    W_A  = diag(diag(covd_A).^-1);
else
    W_ts = diag(ones(N,1));
    W_T = diag(ones(N,1));
    W_A = diag(ones(N,1));
end

% add constraint
G(N+1,:) = 1;
dtstar_pairwise(N+1,:)=0;
dT_pairwise(N+1,:)=0;
dlgA0_pairwise(N+1,:)=0;
W_ts = [[W_ts,zeros(N,1)];[zeros(1,N),1]]; % add a 1 on the end of the diagonal
W_T = [[W_T,zeros(N,1)];[zeros(1,N),1]]; % add a 1 on the end of the diagonal
W_A = [[W_A,zeros(N,1)];[zeros(1,N),1]]; % add a 1 on the end of the diagonal

% kill bad pairs
isbd = (1./misfitnormed_pairwise==0);
G(isbd,:) = [];
dtstar_pairwise(isbd) = [];
dT_pairwise(isbd) = [];
dlgA0_pairwise(isbd) = [];
W_ts(isbd,:) = []; W_ts(:,isbd) = [];
W_T(isbd,:) = []; W_T(:,isbd) = [];
W_A(isbd,:) = []; W_A(:,isbd) = [];

% add a litle damping
eps = 1e-10;
% solve damped weighted LSQR 
delta_tstar = (G'*W_ts*G + eps*eye(Nstas))\G'*W_ts*dtstar_pairwise;
delta_T = (G'*W_T*G + eps*eye(Nstas))\G'*W_T*dT_pairwise;
delta_lgA0 = (G'*W_A*G + eps*eye(Nstas))\G'*W_A*dlgA0_pairwise;
delta_A0 = exp(delta_lgA0); % unlog

B = (G'*W_ts*G + eps*eye(Nstas))\G'*W_ts;
cov_dtstar = B * diag(diag(W_ts).^-1) * B' ;
std_dtstar = diag(cov_dtstar).^0.5;



end