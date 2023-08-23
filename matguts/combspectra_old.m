function [delta_tstar,delta_T,std_dtstar,pairwise] = combspectra(dat,samprate,parms,ifplot)

%% establish some parms
if nargin < 4 
    ifplot = false;
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

amp2phiwt = parms.inv.amp2phiwt;
fmin      = parms.inv.fmin;
fmax      = parms.inv.fmax;
corr_c_skip = parms.inv.corr_c_skip;
ifwt      = parms.inv.ifwt;

fnq = samprate/2;
dt = 1./samprate;
npt = size(dat,1);

Nstas = size(dat,2);
if Nstas == 1, error('Only 1 station''s data!'), end
N = handshake(Nstas);

%% set up results structures
dtstar_pairwise = zeros(N,1);
dT_pairwise = zeros(N,1);
misfitnormed_pairwise = zeros(N,1);

As_pairwise = zeros(N,Nwds);
phis_pairwise = zeros(N,Nwds);
wts_pairwise = zeros(N,Nwds);
doinds_pairwise = zeros(N,Nwds);

%% prepare filter + cleaning parms
% Make set of period windows for bandpass filter
Tmids = logspace(log10(Tmin),log10(Tmax),Nwds)';
if strcmp(Tw_opt,'scale')
    Twdhs = 0.5*diff(logspace(log10(Tmin/2),log10(2*Tmax),Nwds+1)');
else
    Twdhs = 0.5*Tw_opt(:).*ones(size(Tmids));
end

% Twdhs = 0.5*Tmids;
fmids = 1./Tmids;

flos = 1./(Tmids + Twdhs);
fhis = 1./(Tmids - Twdhs);

%% plot filters
if ifplot
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

%% set up cleaning/taper/window parms
nwin=round((postx+prex)/dt); % window length in samples

wdo1 = tukeywin(npt,2*taperx);

n1=round((pretime-prex)/dt); % first sample in window
wdo2=[zeros(n1,1);tukeywin(nwin,2*taperx); zeros(npt-n1-nwin,1)]; % taperx% tukey window
ibds=[max(1,n1-floor(nwin*2*taperx)), min(npt,n1+nwin+floor(nwin*2*taperx))]; % extend so traces are extra 20% of window on either side
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

hw = waitbar(count/N,'Progress through station-station dA,dphi');

for is1 = 1:Nstas
for is2 = is1+1:Nstas
    waitbar(count/N,hw)
    count = count+1;
    % make elements of eventual G matrix
    ui(2*(count-1)+[1 2]) = count;
    uj(2*(count-1)+[1 2]) = [is1 is2];
    u (2*(count-1)+[1 2]) = [-1 1]; % delta is value of 2 - value of 1
    
    %% RUN THROUGH THE COMB - option to correct cycle skip
    [ As,phis,wts ] = run_comb( dat(:,is1),dat(:,is2),fltinfo,wdo1,wdo2,jbds,dt,pretime,maxphi,corr_c_skip,fmin,ifplot );
    if ~ifwt, wts(iw) = 1;end

    %% QC
    inds = find(fmids<=fmax & abs(phis)<maxphi & sqrt(wts)>minacor);

    if length(inds)<4 % only do if there are at least 4 ok measurements!
        As_pairwise(count,:) = As;
        phis_pairwise(count,:) = phis;
        wts_pairwise(count,:) = wts;
        
        dtstar_pairwise(count) = 0;
        dT_pairwise(count) = 0;
        misfitnormed_pairwise(count) = Inf;
        continue % skip to next station pair
    end
    
    %% calculate dtstar and dT
    % estimate from simultaneous inversion of amp and phase data
%     [ dtstar,dT,A0,misfit,E ] = invert_Aphi_4_dtdtstar( As(inds),phis(inds),fmids(inds),fnq, wts(inds),amp2phiwt);
    res = 100*ones(handshake(4),1);
    parms = 3;
    if ifplot
    figure(69), clf, fprintf('\n')
    end
    for iffs = sum(fmids(inds)<fmin):length(inds)
        doinds = inds(end-iffs+1:end);
        [ dtstar_i,dT_i,A0_i,misfit_i,E_i ] = invert_Aphi_4_dtdtstar( As(doinds),phis(doinds),fmids(doinds), wts(doinds),amp2phiwt);
        res_i = sqrt([amp2phiwt*wts(doinds);wts(doinds)]).*E_i;
        P = 1 - ftest(res,parms,res_i,parms);
        if ifplot
        subplot(211),hold on, plot(iffs,rms(res_i),'o')
        subplot(212),hold on, plot(iffs,P,'or')
        end
        if P > 0 % using a VERY loose criterion for including more data, that the P(improvement) is not negative!
            dtstar = dtstar_i;
            dT = dT_i;
            A0 = A0_i;
            misfit = misfit_i;
            res = E_i;
            if ifplot
            fprintf('Using first %.0f inds\n',iffs)
            end
            dodo = doinds;
        end
    end
    
    %% plots
    if ifplot
        figure(3), clf, set(gcf,'pos',[600 600 500,700])

        subplot(211), hold on
        scatter(fmids,log(As),110*wts,'or','MarkerFaceColor','r')
        scatter(fmids(inds),log(As(inds)),110*wts(inds),'ok','MarkerFaceColor','r')
        scatter(fmids(dodo),log(As(dodo)),100*wts(dodo),'o','MarkerFaceColor','b')
        plot(fmids,log(A0) - pi*fmids*dtstar,'g','Linewidth',1.5)
        plot(fmax*[1 1],[-1 1],'--b')
        xlabel('freq','FontSize',18), ylabel('log(Amp)','FontSize',18)
        title(sprintf('Station %.0f vs. station %.0f',is1,is2),'FontSize',20)
        set(gca,'Xscale','linear','xlim',[0.03 1])

        subplot(212), hold on
        scatter(fmids,phis,110*wts,'or','MarkerFaceColor','r')
        scatter(fmids(inds),phis(inds),110*wts(inds),'ok','MarkerFaceColor','r')
        scatter(fmids(dodo),phis(dodo),100*wts(dodo),'o','MarkerFaceColor','b')
        plot(fmids,(log(fnq) - log(fmids))*dtstar./pi + dT,'g','Linewidth',1.5)
        plot(fmax*[1 1],[-1 1],'--b')  
        set(gca,'Xscale','log','xlim',[0.03 1.],'ylim',maxphi*[-1 1])
        title(sprintf('Misfit/N = %.3f  tstar=%.3f  dT=%.3f',misfit./length(inds),dtstar,dT),'FontSize',20)
        xlabel('log(freq)','FontSize',18), ylabel('phase-shift (s)','FontSize',18)
    end
    
    As_pairwise(count,:) = As;
    phis_pairwise(count,:) = phis;
    wts_pairwise(count,:) = wts;
    doinds_pairwise(count,doinds) = true;

    dtstar_pairwise(count) = dtstar;
    dT_pairwise(count) = dT;
    misfitnormed_pairwise(count) = misfit./length(inds); % weight  will be  1./misfit, normalised by number of datapoints

end
end
delete(hw)

pairwise = struct('dtstar',dtstar_pairwise,'dT',dT_pairwise,'misfit_normed',misfitnormed_pairwise,...
                  'As',As_pairwise,'phis',phis_pairwise,'wts',wts_pairwise,'doinds',doinds_pairwise);

%% solve the least squares problem
G = sparse(ui,uj,u,N,Nstas,2*N);

if ifwt
    W = 1./misfitnormed_pairwise; 
else
    W = ones(N,1);
end

% add constraint
G(N+1,:) = 1;
dtstar_pairwise(N+1,:)=0;
dT_pairwise(N+1,:)=0;
W = diag([W;1]);

% kill bad pairs
isbd = (1./misfitnormed_pairwise==0);
G(isbd,:) = [];
dtstar_pairwise(isbd) = [];
dT_pairwise(isbd) = [];
W(isbd,:) = []; W(:,isbd) = [];

% results
delta_tstar = (G'*W*G)\G'*W*dtstar_pairwise;
delta_T = (G'*W*G)\G'*W*dT_pairwise;

cov_dtstar = ( (G'*G)\G' ) * diag(diag(W).^-2) * ( (G'*G)\G' )' ;
std_dtstar = diag(cov_dtstar).^0.5;


end