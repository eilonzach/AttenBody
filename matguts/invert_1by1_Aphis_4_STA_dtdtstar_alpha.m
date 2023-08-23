function [ dtstar_pref,dT_pref,A0_pref,alpha_pref,alpha_misfits,dtstars,dTs,A0s,misfits_amp,misfits_phi ] ...
    = invert_1by1_Aphis_4_STA_dtdtstar_alpha  ( Amat,phimat,fmids,test_alphas,wtmat,amp2phiwt )
%[ dtstar_pref,dT_pref,A0_pref,alpha_pref,alpha_misfits,dtstars,dTs,A0s, misfits_amp, misfits_phi ] ...
%     = allin1invert_Aphis_4_STA_dtdtstar_alpha ( Amat,phimat,freqs,test_alphas,wtmat,amp2phiwt )
%   Script to simultaneously invert pairwise frequency and phase spectra
%   for dtstar and dT at a whole array of stations, looping over a range of
%   alpha values, solving for the best-fitting alpha and the corresponding
%   dtstar and dT
% 
%  Amat and phimat are matrices each of size [handshake(Nstas) x Nfreq]
%  fmids is a vector of frequencies
%  test_alphas is the vector of alpha values (?0) to test
%  wtmat is a matrix of weights, the same size as Amat
% 
% if Q is frequency independent (alpha==0)
%     A = A0*exp(-pi*dtstar*f)
%     ln(A) = ln(A0) - pi*dtstar*f 
%     phi = (ln(f) - ln(fNq))*dtstar/pi + dT
% 
% elseif Q is frequency dependent (alpha~   =0)
%     A = A0*exp(-pi * f0^alpha * f^(1-alpha) * dtstar)
%     ln(A) = ln(A0) - pi * f0^alpha * f^(1-alpha) * dtstar0
%     phi = 0.5*cot(alpha*pi/2)*(f/f0)^alpha + dT


%% prelims
if nargin < 5 || isempty(wtmat)
    wtmat = ones(size(Amat));
end
if nargin < 6 || isempty(amp2phiwt)
    amp2phiwt = 1;
end

Npair = size(Amat,1);
Nstas = quadratic_solve(1,-1,-2*Npair);
Na = length(test_alphas);

%% results structures

dtstars = zeros(Nstas,Na);
dTs = zeros(Nstas,Na);
A0s = zeros(Nstas,Na);
misfits_amp = zeros(Npair,Na);
misfits_phi = zeros(Npair,Na);

% loop on alphas
for ia = 1:length(test_alphas)
    alp = test_alphas(ia);


    % loop over station pairs
    ui = zeros(2*Npair,1);
    uj = zeros(2*Npair,1);
    u  = zeros(2*Npair,1);
    count = 0;

    dtstar_pairwise = zeros(Npair,1);
    dT_pairwise = zeros(Npair,1);
    lgA0_pairwise = zeros(Npair,1);
    misfitnormed_pairwise = zeros(Npair,1);

    for is1 = 1:Nstas
    for is2 = is1+1:Nstas
        count = count+1; % count is the same as the handshake #
        % make elements of eventual G matrix
        ui(2*(count-1)+[1 2]) = count;
        uj(2*(count-1)+[1 2]) = [is1 is2];
        u (2*(count-1)+[1 2]) = [-1 1]; % delta is value of 2 - value of 1

        [ dtstar,dT,A0,misfit,~,m_amp,m_phi ] ...
            = invert_1pair_Aphi_4_dtdtstar(Amat(count,:),phimat(count,:),fmids, wtmat(count,:),amp2phiwt,alp);

        dtstar_pairwise(count) = dtstar;
        dT_pairwise(count) = dT;
        lgA0_pairwise(count) = log(A0);
        misfits_amp(count,ia) = m_amp;
        misfits_phi(count,ia) = m_phi;
        misfitnormed_pairwise(count) = misfit./sum(wtmat(count,:)~=0); % weight  will be  1./misfit, normalised by number of datapoints
    end 
    end

    %% solve the least squares problem
    G = sparse(ui,uj,u,Npair,Nstas,2*Npair);

    W = 1./misfitnormed_pairwise; 

    % add constraint
    G(Npair+1,:) = 1;
    dtstar_pairwise(Npair+1,:)=0;
    dT_pairwise(Npair+1,:)=0;
    lgA0_pairwise(Npair+1,:)=0;
    W = diag([W;1]);

    % kill bad pairs
    isbd = (1./misfitnormed_pairwise==0);
    G(isbd,:) = [];
    dtstar_pairwise(isbd) = [];
    dT_pairwise(isbd) = [];
    lgA0_pairwise(isbd) = [];
    W(isbd,:) = []; W(:,isbd) = [];

    % results
    dtstars(:,ia) = (G'*W*G)\G'*W*dtstar_pairwise;
    dTs(:,ia) = (G'*W*G)\G'*W*dT_pairwise;
    A0s(:,ia) = exp((G'*W*G)\G'*W*lgA0_pairwise);
    

end %loop on alphas

f = amp2phiwt/(amp2phiwt+1);
alpha_misfits = f*sum(misfits_amp,1) + (1-f)*sum(misfits_phi,1);
alpha_misfits = alpha_misfits(:);



%% minimise misfit
alpha_pref = test_alphas(mindex(alpha_misfits));
dtstar_pref = dtstars(:,mindex(alpha_misfits));
dT_pref = dTs(:,mindex(alpha_misfits));
A0_pref = A0s(:,mindex(alpha_misfits));

% figure(77), clf;
% plot(test_alphas,alpha_misfits,'-o')
% xlabel('F-dependency ($\alpha$)','interpreter','latex','FontSize',22)
% ylabel('Global misfit, ($\chi^2$)','interpreter','latex','FontSize',22)
% set(gca,'FontSize',14,'box','on')
% 


end

