function [ Amat,phimat,dtst_pair,dT_pair ] = pred_Amat_phimat( dtstars,dTs,A0s,fmids,alp )
% [ Amat,phimat,dtst_pair,dT_pair ] = pred_Amat_phimat( dtstars,dTs,A0s,fmids,[alpha=0] )
%   function to compute pairwise amplitude and phase spectra, given known
%   values of dtstar, dT, A0 for each station. 
if nargin < 5 || isempty(alp)
    alp = 0;
end

Nstas = length(dtstars);
Npair = handshake(Nstas);
Nf = length(fmids);

Amat = zeros(Npair,Nf);
phimat = zeros(Npair,Nf);
dtst_pair = zeros(Npair);
dT_pair = zeros(Npair);

ws = 2*pi*fmids;

f0=1;
w0 = 2*pi*f0;

count = 0;
for is1 = 1:Nstas
for is2 = is1+1:Nstas
    count = count+1;
    
    RA0 = A0s(is2)/A0s(is1);
    dtst = dtstars(is2)-dtstars(is1);
    dT = dTs(is2)-dTs(is1);
    
    if alp == 0;
        lnApred   = log(RA0) - pi*fmids*dtst;
        phipred = -(1/pi)*log(fmids/f0)*dtst + dtst + dT;
    else
        lnApred = log(RA0) - 0.5*ws.^(1-alp) * w0^alp * dtst;
        phipred = 0.5*cot(alp*pi/2) * (ws/w0).^(-alp) * dtst + dT;
    end
    
    Amat(count,:) = exp(lnApred);
    phimat(count,:) = phipred;    
    dtst_pair(count,:) = dtst;
    dT_pair(count,:) = dT;


end
end



end

