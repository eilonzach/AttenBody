function [delta_tstar,delta_T,std_dtstar,pairwise,ff_use,delta_A0] = fAxPspectra(dat,fcross,samprate,parms,ifplot,mtmspecss,slatlon)
% [delta_tstar,delta_T,std_dtstar,pairwise,fmids,delta_A0] = fAxPspectra(dat,fcross,samprate,parms,ifplot,mtmspecss,slatlon)

%% establish some parms
if nargin < 5 
    ifplot = false;
end
if nargin < 6
    parms.qc.mtmRcomp = 0;
    mtmspecss = [];
end
if nargin <7
    parms.qc.maxdist = false;
    slatlon = [];
end

dt = 1./samprate;

pretime = parms.wind.pretime;
prex    = parms.wind.prex;
postx   = parms.wind.postx;
taperx  = parms.wind.taperx;
fhi     = parms.wind.fhi;
flo     = parms.wind.flo;

% do a little calculation to account for removal of taper on prex and postx
% we want to keep everything truly within prex and postx, so add a little
% time on either side for the taper
padt = taperx*(postx+prex)./(1 - 2*taperx);
postx = postx + padt;
prex = prex + padt;

minacor = parms.qc.minacor;
maxphi  = parms.qc.maxphi;
try
    maxdA = parms.qc.maxdA;
catch
    maxdA = 0.3;
end

amp2phiwt = parms.inv.amp2phiwt;
ifunwrap = parms.fAxP.ifunwrap;
% fmax      = parms.inv.fmax;

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

if ~isfield(parms.qc,'mtmRcomp') % specify whether to compare amplitudes to specR
    parms.qc.mtmRcomp = false;
elseif isempty(parms.qc.mtmRcomp)
    parms.qc.mtmRcomp = false;
end

if ~isfield(parms.qc,'phi2deriv') % specify whether to use phi second derivative to weight spectra
    parms.qc.phi2deriv = false;
elseif isempty(parms.qc.phi2deriv)
    parms.qc.phi2deriv = false;
end

if ~isfield(parms.qc,'maxdist') % specify whether to compare amplitudes to specR
    parms.qc,maxdist = false;
elseif isempty(parms.qc.maxdist)
    parms.qc,maxdist = false;
end



%% filter and window all the traces
% do NOT want to norm here
cp = struct('samprate',1./dt,'pretime',pretime,'prex',prex,'postx',postx,...
                'taperx',taperx,'fhi',fhi,'flo',flo,'npoles',2,'norm',0,'detrend',1);
[ datwf ] = data_clean( dat,cp); 

% clean the data (shouldn't need it!)
% nnan=~isnan(datwf);
% datwf = reshape(detrend(datwf(nnan)),size(datwf)); % detrend non-nan data
% datwf(isnan(datwf)) = 0; %set any nans to zero


%% Calculate the complex fourier spectra for all stations
[datwf_fft,ff_fft] = fft_ze(datwf,dt);
% make one-sided, non-inclusive of zero
f_ipos = (ff_fft>0);
datwf_fft = datwf_fft(f_ipos,:);
ff_fft = ff_fft(f_ipos);

% subset to widest frequency range we can use
f_iuse = find((ff_fft>=min(fcross(:,1))) & (ff_fft<=max(fcross(:,2))));

ff_use = ff_fft(f_iuse);

% whiletesting compare spectra from MTM and from fft
% figure(443); clf, hold on
% plot(ff_fft,abs(datwf_fft),'-','linewidth',2)
% plot(mtmspecss.frq,mtmspecss.specs,'--','linewidth',1.5)
% set(gca,'xscale','log')

%% prep for results
fnq = samprate/2;
Nf = length(f_iuse);

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

% better to make these things sparse
As_pairwise = sparse(N,Nf);
phis_pairwise = sparse(N,Nf);
wts_pairwise = sparse(N,Nf);
inds_pairwise = sparse(N,Nf);



% N = length(dat1_fft);
% fplt = (ff_fft>=max(fcross(:,1))) & (ff_fft<=min(fcross(:,2)));
% specRfft = smooth(log(abs(dat2_fft)./abs(dat1_fft)));
% % specPfft = smooth(unwrap(angle(qdatwc2_fft./dat1_fft))./(-2*pi*ff_fft));
% specPfft = smooth((angle(dat2_fft./dat1_fft))./(-2*pi*ff_fft));
% 
% specRfft = specRfft(fplt);
% specPfft = specPfft(fplt);
% ff_fft = ff_fft(fplt);
% % wt_fft = 1./(0.5+abs(interp1(eqar(ind1).frq,log(eqar(ind2).specs./eqar(ind1).specs),ff_fft) - specRfft));
% wt_fft_a = 0.4./(.2+abs(interp1(eqar(ind1).frq,log(eqar(ind2).specs./eqar(ind1).specs),ff_fft) - specRfft));
% wt_fft_p = [1;abs(diff(specPfft)./diff(ff_fft))<10];
% wt_fft = (wt_fft_a.*wt_fft_p) + 0.0001;
% figure(444),clf
% subplot(211), hold on
% plot(ff_fft,specRfft ,'-b','linewidth',2)
% plot(eqar(ind1).frq,log(eqar(ind2).specs./eqar(ind1).specs),'-r')
% plot(ff_fft,wt_fft,'k')
% xlim([0 0.5])
% subplot(212), hold on
% plot(ff_fft,specPfft ,'-b','linewidth',2)
% xlim([0 0.5])



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

    %% skip immediately if stations too far apart
    if parms.qc.maxdist
        interstadist = distance(slatlon(is1,:),slatlon(is2,:));
        if interstadist>parms.qc.maxdist % only do stations closer than critical distance
%             As_pairwise(count,:) = sparse(1,Nf);
%             phis_pairwise(count,:) = sparse(1,Nf);
%             wts_pairwise(count,:) = sparse(1,Nf);
            
            dtstar_pairwise(count) = nan;
            dT_pairwise(count) = nan;
            misfitnormed_pairwise(count) = Inf;
            continue % skip to next station pair
        end
    end

    DAT1 = datwf_fft(:,is1); % use the full f-range of the spectra
    DAT2 = datwf_fft(:,is2);

    %% Calculate the differential spectra
    As = exp(smooth(log(abs(DAT2)./abs(DAT1)))); % NOT in log space, but smoothed in log space
    if ifunwrap
        phis = smooth(unwrap(angle(DAT2./DAT1))./(-2*pi*ff_fft));
    else
        phis = smooth((angle(DAT2./DAT1))./(-2*pi*ff_fft));
    end 
    %  subset to usable frequencies AFTER the smoothing
    As = As(f_iuse);
    phis = phis(f_iuse);

    %% calculate the weightings
    % wt_fft = 1./(0.5+abs(interp1(eqar(ind1).frq,log(eqar(ind2).specs./eqar(ind1).specs),ff_fft) - specRfft));
    
    % opt to compare the measured amplitude spectrum to the MTM 
    if parms.qc.mtmRcomp
        As_mtmR = interp1(mtmspecss.frq,mtmspecss.specs(:,is2)./mtmspecss.specs(:,is1),ff_use);
        wts_mtmRcomp = .5./(0.25 + abs(log(As) - log(As_mtmR))); % water level it to 0.5. 
        wts_mtmRcomp = wts_mtmRcomp(:);
        % if there is ZERO difference in log amplitudes, then gets double
        % weighted. If there is 0.25 ln unit difference, gets normal weight.
        % If there is 1.25 difference in weight, down to third weight
    else
        wts_mtmRcomp = ones(size(ff_use));
    end

%         % while testing
%         figure(444),clf
%         subplot(211), hold on
%         plot(ff_use,log(As) ,'-ob','linewidth',2)
%         plot(ff_use,log(As_mtmR),'-or')
%         % plot(ff_fft,wt_fft,'k')
%         xlim([0 0.5])
%         subplot(212), hold on
%         plot(ff_use,phis ,'-ob','linewidth',2)
%         xlim([0 0.5])
%         pause
%         end
%         end
%         % end while testing

    % weight the phase spectrum by whether or not the second derivatives are stupid
    % if low gradient, full weight. As gradient get above, downweight
    if parms.qc.phi2deriv
        obs2ndderiv = gradient(gradient(phis,ff_use),ff_use); % observed second derivative
        max2ndderiv = abs(gradient(gradient((3.5/pi)*log(ff_use),ff_use),ff_use)); % theorietical second derivative for very high dtstar
        wts_phi2deriv = 4*max2ndderiv./(2*max2ndderiv + abs(max2ndderiv - abs(obs2ndderiv))); % water level it to max2ndderiv. 
        wts_phi2deriv = sqrt(wts_phi2deriv(:));

    else
        wts_phi2deriv = ones(size(ff_use));
    end
    wts = (wts_mtmRcomp.*wts_phi2deriv);
    
    
    if ~ifwt, wts(:) = 1; end
    wts(As<0) = 0;
    As(As<0) = 0;
    wts(phis==maxphi) = 0;
    
    % work out frequency limit for this pair
    fmin_pair = max(fcross([is1,is2],1));
    fmax_pair = min(fcross([is1,is2],2));


    %% QC
    inds = (ff_use>=fmin_pair) & ...
           (ff_use<=fmax_pair) & ...
           (abs(phis)<maxphi);
       
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
    
	% re-norm wts temporarily
    wtnormfac = rms(wts(inds));
    wts(inds) = wts(inds)./wtnormfac; % DANGER - if all wts zero, gives inf. Should not have all zero - see above if statement.
    if any(isinf(wts))
        keyboard
    end
    if any(isinf(phis)) || any(isinf(As)) 
        keyboard
    end

    
    %% calculate dtstar and dT
    % estimate from simultaneous inversion of amp and phase data
    [ dtstar,dT,A0,misfit,resnorm,A_SSE,P_SSE,lnA_SSE,dtstar_std,dT_std,dlnA0_std ] ...
        = invert_1pair_Aphi_4_dtdtstar( As(inds),phis(inds),ff_use(inds), wts(inds),amp2phiwt,alpha);
    
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

    % un-normalise weights (for comparison among station pairs)
    wts(inds) = wts(inds).*wtnormfac;
    
    %% plots
    if ifplot~=0
        figure(4), clf, set(gcf,'pos',[600 600 500,700])

        subplot(211), hold on
        pwts(pwts<=0)=nan;
        scatter(ff_use,log(As),110*pwts,'or','MarkerFaceColor','r')
        scatter(ff_use(inds),log(As(inds)),110*pwts(inds),'ok','MarkerFaceColor','r')
%         scatter(fmids(dodo),log(As(dodo)),100*wts(dodo),'o','MarkerFaceColor','b')
        if alpha==0
            plot(ff_use,log(A0) - pi*ff_use*dtstar,'g','Linewidth',1.5)
        else
            plot(ff_use,log(A0) - pi*ff_use.^(1-alpha)*dtstar,'g','Linewidth',1.5)
        end
        plot(fmin_pair*[1 1],[-1 1],'--b')
        plot(fmax_pair*[1 1],[-1 1],'--b')
        xlabel('freq','FontSize',18), ylabel('log(Amp)','FontSize',18)
        title(sprintf('Station %.0f vs. station %.0f',is1,is2),'FontSize',20)
        set(gca,'Xscale','linear','xlim',[0.03 1])
        if ~isempty(mtmspecss)
            plot(mtmspecss.frq,log(mtmspecss.specs(:,is2)./mtmspecss.specs(:,is1)),'-xk')
        end

        subplot(212), hold on
        scatter(ff_use,phis,110*pwts,'or','MarkerFaceColor','r')
        scatter(ff_use(inds),phis(inds),110*pwts(inds),'ok','MarkerFaceColor','r')
%         scatter(fmids(dodo),phis(dodo),100*wts(dodo),'o','MarkerFaceColor','b')
        if alpha==0
            plot(ff_use,(log(fnq) - log(ff_use))*dtstar./pi + dT,'g','Linewidth',1.5)
        else
            plot(ff_use,0.5 * cot(alpha*pi/2) * ff_use.^-alpha * dtstar + dT,'g','Linewidth',1.5)
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