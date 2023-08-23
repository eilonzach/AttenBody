function [dcor,dcstd,dvcstd,acor]=xcortimes(dtrn,dt,pretime,lagmax,iplot)
% [dcor, dcstd, dvcstd, acor]=xcortimes(dtrn, dt, pretime, lagmax,iplot)
%         something sort of like the vandecar and crosson cross correlation
%  assumes have array dtrn = (npt x nsta) 
%  dt = samplerate
%  pretime = time of first sample relative to 0 lag
%  lagmax = max allowed lag in seconds to test
%  iplot=1 to do graphics, 0 to just return lags
%
% dcstd = Geoff's version of the stds
% dvcstd = stds from Vandecar & Crosson eqn 8
% acor = mean correlation coefficient of each trace (with stack)
% 
% NEGATIVE DCOR MEANS EARLIER ARRIVAL, POSITIVE MEANS LATER 
% 
% edits ZE 08/2013

if nargin < 3 || isempty(pretime)
    pretime = size(dtrn,1)*dt/2;
end
if nargin < 4 || isempty(lagmax)
    lagmax = size(dtrn,1)*dt/10;
end
if nargin < 5 || isempty(iplot)
    iplot = 0;
end


nsta=size(dtrn,2);
npt=size(dtrn,1);
nmat=nsta*(nsta-1)/2; % total number of combinations of statinos
Gmat=zeros(nmat+1,nsta);
avec=zeros(nmat+1,1);
lmax=round(lagmax./dt);   % Max lag allowed, should set in fcn
kk=0;
%if (iplot==1) disp('building lag matrix...'); end
for ii=1:(nsta-1) % for each station...
    for jj=(ii+1):nsta % loop over combos with all stations after it
        kk=kk+1;
        [c,lags]=xcorr(dtrn(:,ii),dtrn(:,jj),lmax,'unbiased');
        pk=find(c==max(c));
        avec(kk) = lags(pk(1));
        Gmat(kk,ii) = 1;
        Gmat(kk,jj) = -1;
    end
end
Gmat(kk+1,:)=ones(1,nsta); % damping
%if (iplot==1) disp('Inverting lags...'); end
dcor = Gmat \ avec;
resid=avec - Gmat*dcor;
dcstd=sqrt(diag(inv(Gmat'*Gmat)).*var(resid));

% Find standard dev of residuals using Vandecar & Crosson equation 8
R = zeros(nsta) ;
R(tril(true(nsta),-1)) = resid(1:end-1);
R = R - R';
V = R.^2; 
dvcstd = sqrt(sum(V,2)./(nsta-2));

% Find acor
stak = zeros(npt,1);
tstak = [0:npt-1]';
for is=1:nsta
    stak=stak+interp1(tstak-dcor(is),dtrn(:,is),tstak);
end
stak=stak./nsta;
jmx=npt-ceil(max(dcor));
jmn=ceil(abs(min(dcor)))+1;
ida=find(isfinite(stak) & (1:npt)'<jmx & (1:npt)'>jmn );
acor = zeros(nsta,1);
for is=1:nsta
    acor(is)=sum(stak(ida).*dtrn(ida+round(dcor(is)),is))./(std(stak(ida),1).*std(dtrn(ida+round(dcor(is)),is),1))./length(ida);
end

dcor = dcor .* dt;
dcstd=dcstd .* dt;

if (iplot==1) 
    figure(46)
    clf
    %disp('...finishing up');
    tt0=dt.*(0:(npt-1))-pretime; %edited zje: was tt0=dt.*(0:(npt-1))-pretime;
for ii=1:nsta
    subplot(211);
    plot(1:npt,dtrn(:,ii),'Linewidth',1.5);
    title('Pre-alignment')
    hold on;
    subplot(212);
    plot(tt0-dcor(ii),dtrn(:,ii),'Linewidth',1.5);
    title('Post-alignment')
    hold on;
end
end
return
