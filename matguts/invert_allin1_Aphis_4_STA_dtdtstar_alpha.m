function [ dtstar_pref,dT_pref,A0_pref,alpha_pref,alpha_misfits,...
             dtstars,dTs,A0s,misfits_amp,misfits_phi,dtstars_std] ...
            = invert_allin1_Aphis_4_STA_dtdtstar_alpha( Amat,phimat,fmids,test_alphas,wtmat,amp2phiwt,w0 )
%[ dtstar_pref,dT_pref,A0_pref,alpha_pref,alpha_misfits,
%   dtstars,dTs,A0s, misfits_amp, misfits_phi,dtstars_std] ...
%    = invert_allin1_Aphis_4_STA_dtdtstar_alpha( Amat,phimat,fmids,test_alphas,wtmat,amp2phiwt,w0 )
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
%     A = A0*exp(-(pi/((2*pi)^alpha)) * f^(1-alpha) * dtstar)
%     ln(A) = ln(A0) - (pi/((2*pi)^alpha)) * f^(1-alpha) * dtstar
%     phi = 0.5*cot(alpha*pi/2)*f^alpha + dT


%% prelims
if nargin < 5 || isempty(wtmat)
    wtmat = ones(size(Amat));
end
if nargin < 6 || isempty(amp2phiwt)
    amp2phiwt = 1;
end
if nargin < 7 || isempty(w0)
    w0 = 2*pi;
end

Npair = size(Amat,1);
Nstas = quadratic_solve(1,-1,-2*Npair);
Nf = length(fmids);
Na = length(test_alphas);

ws = 2*pi*fmids;


%% results structures
dtstars = zeros(Nstas,Na);
dTs = zeros(Nstas,Na);
A0s = zeros(Nstas,Na);

dtstars_std = zeros(Nstas,Na);

wts = reshape(wtmat',Npair*Nf,1);

d_Amp = reshape(log(Amat)',Npair*Nf,1);
d_phi = reshape(phimat',Npair*Nf,1);

% seek infinite data
if any(isinf(d_phi))
    warning('Some infinite phi values - zeroing')
    wts(isinf(d_phi)) = 0;
    d_phi(isinf(d_phi)) = 0;
end
if any(isinf(d_Amp))
    warning('Some infinite amp values - zeroing')
    wts(isinf(d_Amp)) = 0;
    d_Amp(isinf(d_Amp)) = 0;
end


G_Amp = spalloc(Npair*Nf,3*Nstas,4*Nf*Npair);
G_phi = spalloc(Npair*Nf,3*Nstas,4*Nf*Npair);

alpha_misfits = zeros(Na,1);

misfits_amp = nan(Npair,Na);
misfits_phi = nan(Npair,Na);

for ia = 1:length(test_alphas)
    alpha = test_alphas(ia);
    
    if alpha==0
        Ax = -0.5*ws;
        Px = 1 - log(ws/2/pi)./pi;
    else
        Ax = -0.5 * w0^alpha * ws.^(1-alpha);
        Px = 0.5 * cot(alpha*pi/2) * (ws/w0).^(-alpha);
    end
    
    % data vector: [Npair*Nf x 1] for both A and phi
    % mod parm vector: [A0s(Nstas); dtstar(Nstas); dT(Nstas)]
    
    if any(isinf(Ax)) || any(isinf(Px))
        keyboard
    end
    
    count = 0;
    for is1 = 1:Nstas
    for is2 = is1+1:Nstas
        count = count+1; % count is the same as the handshake #

        % slot into G matrices
        yind = [1:Nf] + (count-1)*Nf;

        G_Amp(yind,is1)         = -1;
        G_Amp(yind,is2)         = 1;
        G_Amp(yind,Nstas+is1)   = -Ax;
        G_Amp(yind,Nstas+is2)   = Ax;

        G_phi(yind,Nstas+is1)   = -Px;
        G_phi(yind,Nstas+is2)   = Px;
        G_phi(yind,2*Nstas+is1) = -1;
        G_phi(yind,2*Nstas+is2) = 1;
    end 
    end

    constraint_A0  = sparse(1,[1:Nstas]        ,1,1,3*Nstas);
    constraint_Amp = sparse(1,[1:Nstas]+Nstas  ,1,1,3*Nstas);
    constraint_phi = sparse(1,[1:Nstas]+2*Nstas,1,1,3*Nstas);
    
    d_all = [d_Amp;d_phi;0;0;0];
    G_all = [G_Amp;G_phi;constraint_A0;constraint_Amp;constraint_phi];
    
    w_all = [amp2phiwt*wts;wts;1;1;1];
    N = length(w_all);
    spdw_all = sparse(1:N,1:N,w_all,N,N,N);

    M = size(G_all,2);
    
    dw = spdw_all*d_all;
    Gw = spdw_all*G_all;
    
    % kill allzero rows
%     zerorows = sum(abs(Gw),2)==0;
%     dw = dw(~zerorows);
%     Gw = Gw(~zerorows,:);
    
    %% try to solve uniquely
%     m = (Gw'*Gw)\Gw'*dw;   
    try
        m = lsqr(Gw,dw,1e-5,200);
        if isnan(m), error('lsqr failed. trying some small damping'); end
    catch
        M = 3*Nstas;
        H = speye(M)*1e3;
        h = sparse(M,1);
        F = [Gw;H];
        f = [dw;h];
        m = (F'*F)\F'*f
%         m = lsqr(F,f,1e-1,500);
    end
    
    
    dtstars(:,ia) = m(Nstas+ [1:Nstas]);
    dTs(:,ia) = m(2*Nstas + [1:Nstas]);
    A0s(:,ia) = exp(m(1:Nstas));

    E =  [G_all*m - d_all];
    Ew = [Gw*m - dw];
 
    %% estimate model error from posteriori data error
    dof = Npair*Nf - M;
    iamp = [1:Npair*Nf]';
    iphi = [1:Npair*Nf]' + Npair*Nf;
    % sd_posterior = sqrt(sum(e^2)/dof) % like RMS but normalised by d.o.f.
    sd_amp_posteriori = sqrt(sum(E(iamp).^2)./dof);
    sd_phi_posteriori = sqrt(sum(E(iphi).^2)./dof);
    % weighted average (weight already in numerator in Ew, add to denom
    dof_wamp = full(sum(diag(spdw_all(iamp,iamp).^2)) - M*mean(diag(spdw_all(iamp,iamp).^2)));
    dof_wphi = full(sum(diag(spdw_all(iphi,iphi).^2)) - M*mean(diag(spdw_all(iphi,iphi).^2)));
    sdw_amp_posteriori = sqrt(sum(Ew(iamp).^2)./dof_wamp);
    sdw_phi_posteriori = sqrt(sum(Ew(iphi).^2)./dof_wphi);



    covd = spdiags([sdw_amp_posteriori.^2*ones(Npair*Nf,1);...
                  sdw_phi_posteriori.^2*ones(Npair*Nf,1);...  
                  0;0;0],0,2*Npair*Nf+3,2*Npair*Nf+3);
    Gg = (Gw'*Gw)\Gw'*spdw_all;
%     R = (Gw'*Gw)\(Gw'*Gw); 
%     N = Gw*Gg;
    covm = Gg*covd*Gg';
    varm = full(diag(covm));

    dtstars_std(:,ia) = sqrt(varm(Nstas+ [1:Nstas]));
    dTs_sd = sqrt(varm(2*Nstas+ [1:Nstas])) ;   
    A0s_sd = sqrt(varm([1:Nstas]));

    %% misfits
    alpha_misfits(ia) = E'*spdw_all*E;
    
    misfits_amp(:,ia) = sum(reshape(wts.*E(           [1:Npair*Nf]).^2,Nf,Npair))'./sum(reshape(wts,Nf,Npair))';
    misfits_phi(:,ia) = sum(reshape(wts.*E(Npair*Nf + [1:Npair*Nf]).^2,Nf,Npair))'./sum(reshape(wts,Nf,Npair))';
end

%% minimise misfit
alpha_pref = test_alphas(mindex(alpha_misfits));
dtstar_pref = dtstars(:,mindex(alpha_misfits));
dT_pref = dTs(:,mindex(alpha_misfits));
A0_pref = A0s(:,mindex(alpha_misfits));

if length(test_alphas)>1
figure(77), clf;
plot(test_alphas,alpha_misfits,'-o')
xlabel('F-dependency ($\alpha$)','interpreter','latex','FontSize',22)
ylabel('Global misfit, ($\chi^2$)','interpreter','latex','FontSize',22)
set(gca,'FontSize',14,'box','on')
end

end

