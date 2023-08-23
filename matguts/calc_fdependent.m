function [ delta_tstar_pref,delta_T_pref,A0_pref,alpha_pref,...
           alpha_VR,alpha_Ew,misfitsw_amp,misfitsw_phi,delta_tstar_std_pref ] ...
    = calc_fdependent( Amat,phimat,fmids,test_alphas,wtmat,amp2phiwt,opt,titlstr,ifplot )
% [ delta_tstar_pref,delta_T_pref,A0_pref,alpha_pref,...
%    alpha_VR,alpha_Ew,misfitsw_amp,misfitsw_phi,delta_tstar_std_pref ]
%    = calc_fdependent( Amat,phimat,fmids,test_alphas,wtmat,[amp2phiwt=2],[opt=1],[titstr=''],[ifplot=0])
%   Function to re-do all the calculations of dtstar and dT testing the
%   various frequency-dependent cases. 
% 
%  Opt==1 ==> all in one
%  Opt==2 ==> one-by-one
if nargin < 6 || isempty(amp2phiwt)
    amp2phiwt = 2;
end
if nargin < 7 || isempty(opt)
    opt = 1;
end
if nargin < 8 || isempty(titlstr)
    titlstr = '';
end
if nargin < 9 || isempty(ifplot)
    ifplot = 0;
end

Na = length(test_alphas);
VRwa = zeros(Na,1);
VRwp = zeros(Na,1);
Eaw = zeros(Na,1);
Epw = zeros(Na,1);

% safety things
wtmat(Amat<=0) = 0;
Amat(Amat<=0) = 1e-9;

if opt==1
 %% all in one method
    [ ~,~,~,~,~,dtstars_a,dTs_a,A0s_a,misfitsw_amp,misfitsw_phi,dtstars_std_a] ...
        = invert_allin1_Aphis_4_STA_dtdtstar_alpha(Amat,phimat,fmids,test_alphas,wtmat,amp2phiwt );
    methstr = 'all in one';
elseif opt == 2
   %% one-by-one method
    [ ~,~,~,~,~,dtstars_a,dTs_a,A0s_a,misfitsw_amp,misfitsw_phi ] ...
        = invert_1by1_Aphis_4_STA_dtdtstar_alpha(Amat,phimat,fmids,test_alphas,wtmat,amp2phiwt );
    methstr = 'one-by-one';
end

% plot
if ~ifplot
    try close(77); end
end

for ia = 1:length(test_alphas)
    % compute A and phi predictions for best fitting values at each alpha
    [ Amat_pred,phimat_pred,dtst_pair,dT_pair ] = pred_Amat_phimat( dtstars_a(:,ia),dTs_a(:,ia),A0s_a(:,ia),fmids,test_alphas(ia) ); 
    
    
    veaw = nanvar((log(Amat_pred(:)) - log(Amat(:))).*wtmat(:)); % variance of the error in Amp
    vaw = nanvar(log(Amat(:)).*wtmat(:));                        % variance of the Amp
    VRwa(ia) = 1 - veaw/vaw;                                          % variance reduction in Amp
    
    vepw = nanvar((phimat_pred(:) - phimat(:)).*wtmat(:));       % variance of the error in phi
    vpw = nanvar(phimat(:).*wtmat(:));                           % variance of the phi
    VRwp(ia) = 1 - vepw/vpw;                                          % variance reduction in phi
    
    Eaw(ia) = sqrt(sum(((log(Amat_pred(:)) - log(Amat(:))).^2).*wtmat(:))./sum((log(Amat(:)).^2).*wtmat(:)));
    Epw(ia) = sqrt(sum(((phimat_pred(:) - phimat(:)).^2).*wtmat(:))./sum((phimat(:).^2).*wtmat(:)));
    
    % plot...
%     figure(44), clf, set(gcf,'pos',[ 440 -139 1430 797])
%     figure(45), clf, set(gcf,'pos',[ 440 -139 1430 797])
%     rind = randsample(size(Amat,1),56);
%     figure(44), for ir = 1:length(rind); subplot(8,7,ir), hold on, scatter(fmids,log(Amat(rind(ir),:)),'or','filled'),plot(fmids,log(Amat_pred(rind(ir),:)),'Linewidth',2), title(num2str(dtst_pair(ir))), end
%     figure(45), for ir = 1:length(rind); subplot(8,7,ir), hold on, scatter(fmids,phimat(rind(ir),:),'og','filled'),plot(fmids,phimat_pred(rind(ir),:),'Linewidth',2), title(num2str(dtst_pair(ir))), end
end
f = amp2phiwt/(amp2phiwt+1);
alpha_Ew = f*Eaw + (1-f)*Epw;
alpha_VR = f*VRwa + (1-f)*VRwp;


if ifplot
figure(79), clf, set(gcf,'pos',[10 10 500 800])
subplot(311)
plot(test_alphas,Eaw,'-or')
title([titlstr,' ',methstr],'FontSize',22,'interpreter','latex')
ylabel('Eaw','FontSize',18,'interpreter','latex')
set(gca,'XtickLabel',[],'xtick',test_alphas,'fontsize',14)

subplot(312)
plot(test_alphas,Epw,'-ob')
ylabel('Epw','FontSize',18,'interpreter','latex')
set(gca,'XtickLabel',[],'xtick',test_alphas,'fontsize',14)

subplot(313)
plot(test_alphas,alpha_Ew,'-ok')
ylabel('Wt misfit','FontSize',18,'interpreter','latex')
xlabel('alphas','FontSize',18,'interpreter','latex')
set(gca,'xtick',test_alphas,'fontsize',14)
end

ia = mindex(-alpha_VR);

delta_tstar_pref = dtstars_a(:,ia);
delta_T_pref = dTs_a(:,ia);
A0_pref = A0s_a(:,ia);
alpha_pref = test_alphas(ia);
delta_tstar_std_pref = dtstars_std_a(:,ia);

end

