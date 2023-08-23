function [ delta_tstar,cov_dtstar,std_dtstar ] = xspecratio( specs,frq,fmax,fmin,wtopt,plotopt )
%[ delta_tstar ] = xspecratio( specs,frq,fmax,wtopt,plotopt )
%   Function to calculate differential t-star by taking pair-wise spectral
%   ratios for all combinations of stations' spectra and then using a
%   least squares approach to solve for best-fitting delta t-star with the
%   constraint that the mean dtstar=0. 
%   N.B. this script assumes that the spectral ratio can be fit by a linear
%   relationship, and the fit is over the frequency interval 0 to fmax
% 
%  INPUTS:
%    specs   - [nfreq x nstas] matrix of all the spectra, in columns
%    frq     - [nfreq x 1] vector of frequencies 
%    fmax    - high frequency end of fitting window (either a single value
%              for all or one per station - min fmax gets used in pair)
%    fmin    - low frequency end of fitting window (either a single value
%              for all or one per station - max fmin gets used in pair)
%    wtopt   - (optional, default=0) to weight least-squares inversion by
%               formal error in the linear fit to the spectral ratio
%    plotopt - (optional, default=0) to plot spectra over plotting interval
% 
%  OUTPUT:
%    delta_tstar - differential tstar
%    cov_dtstar  - estimated covariance matrix for dtstar - standard devs
%                   estimated by std_dstar = diag(cov_dtstar).^-0.5
%    std_dstar   - estimated standard devs of dtstar, as above

if nargin < 4
    wtopt = 0;
end
if nargin < 5
    plotopt = 0;
end

M = size(specs,2);
N = handshake(M);

if length(fmax)==1
    fmax = fmax.*ones(M,1);
end
if length(fmin)==1
    fmin = fmin.*ones(M,1);
end

dtst = zeros(N,1);
wtst = zeros(N,1);

ui = zeros(2*N,1);
uj = zeros(2*N,1);
u  = zeros(2*N,1);
count = 0;

if plotopt
    figure(44); clf; 
end



for ii = 1:M
for jj = ii+1:M
    % spectral ratio
    R = specs(:,ii)./specs(:,jj);
    % take logarithm
    lnR = log(R); 
        
    fcrosslo = max(fmin([ii,jj]));
    fcrosshi = min(fmax([ii,jj]));
    ind = (frq >=  fcrosslo)  &  (frq <= fcrosshi);

    if sum(ind) < 4 % fewer than 3 frequencies usable
        continue
    end

    fo = fit(frq(ind),lnR(ind),'poly1');
    cint = confint(fo)/pi;
    
    % we have a result... add to count
    count = count+1;

    dtst(count) = fo.p1/(-pi);
    wtst(count) = 0.25*diff(cint(:,1)'); % confidence interval spans 4 standard devs. 
    
    if plotopt
    hold on
    hp = plot(frq(ind),lnR(ind),'o-');
    xlim([fcrosslo-0.1 fcrosshi+0.1])
%     plot(fo)
    end
    
    ui(2*(count-1)+[1 2]) = count;
    uj(2*(count-1)+[1 2]) = [ii jj];
    u (2*(count-1)+[1 2]) = [1 -1];
end
end
N = count;
% get rid of excess in vectors
ui(2*N + 1 : end) = [];
uj(2*N + 1 : end) = [];
u(2*N + 1 : end) = [];
dtst(N + 1 : end) = [];
wtst(N + 1 : end) = [];

% build G
G = sparse(ui,uj,u,N,M,2*N);

if wtopt
    W = diag(wtst.^-1); 
else
    W = eye(N);
end

% add constraint
G(N+1,:) = 1;
dtst(N+1,:)=0;
W = diag([diag(W);1]);

delta_tstar = (G'*W*G)\G'*W*dtst;
cov_dtstar = ( (G'*G)\G' ) * diag(diag(W).^-2) * ( (G'*G)\G' )' ;
std_dtstar = diag(cov_dtstar).^0.5;
end

