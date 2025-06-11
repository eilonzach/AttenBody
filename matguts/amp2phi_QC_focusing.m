function [gdofgd,focusscore] = amp2phi_QC_focusing(parms,pairwise,fmids,ifplot,stas)
%
% 
% function to test for any telltale signs of extrinsic scattering or
% focusing and throw out pairwise measurements/stations accordingly so that
% this effect does not bleed into the t-star *intrinsic attenuation*
% measurements

if nargin < 4
    ifplot = false;
end
if nargin < 5
    stas = '';
end

focus_thresh = parms.qc.focus_threshold;
test_a2pweights = [50,5,1]; % will try these in both directions, i.e. 1/these, also

% parse inputs
Amat = pairwise.As;
phimat = pairwise.phis;
wtmat = double(pairwise.inds).*pairwise.wts;

a2pwt_default = parms.inv.amp2phiwt;

% not the place to test alphas...
if isempty(parms.inv.alpha)
    parms.inv.alpha = 0;
end

% primary part of test is the consistency between the tstar estimates from
% amplitude and phase spectra
  test_a2pwt = sort(unique([a2pwt_default,test_a2pweights,1./test_a2pweights]')); %#ok<TRSRT> 

[~,iallamp] = max(test_a2pwt);
[~,iallphi] = min(test_a2pwt);


    Ntest = length(test_a2pwt);
    Npair = size(Amat,1);
    % work out pair# to sta-combination...
    Nstas_dtstar = unhandshake(Npair);
    ipair = 0;
    is1is2 = nan(size(Amat,1),2);
    % work out pairs
    for is1 = 1:Nstas_dtstar
        for is2 = is1+1:Nstas_dtstar
            ipair = ipair+1; 
            is1is2(ipair,:) = [is1,is2]; 
        end
    end

    %% actually do the test, running through different amp2phi wt values
    % prep results
    a2ptest_dtstar = zeros(Nstas_dtstar,Ntest);
    a2ptest_dtstarstd = zeros(Nstas_dtstar,Ntest);
    a2ptest_dT = zeros(Nstas_dtstar,Ntest);
    a2ptest_A0 = zeros(Nstas_dtstar,Ntest);
    a2ptest_VR = zeros(Ntest,1);
    a2ptest_Ew = zeros(Ntest,1);
    a2ptest_Eamp = zeros(Npair,Ntest);
    a2ptest_Ephi = zeros(Npair,Ntest);

    % run through a2pwts
    for itest = 1:length(test_a2pwt)
        fprintf('Test a2pwt %6.2f  ',test_a2pwt(itest))
            [ a2ptest_dtstar(:,itest),a2ptest_dT(:,itest),a2ptest_A0(:,itest),...
                ~,a2ptest_VR(itest),a2ptest_Ew(itest),a2ptest_Eamp(:,itest),a2ptest_Ephi(:,itest),a2ptest_dtstarstd(:,itest) ] ...
            = calc_fdependent( Amat,phimat,fmids,parms.inv.alpha,wtmat,test_a2pwt(itest),parms.inv.opt,[],0);
    end

    % station average errors (add where that station is in a pair)
    a2ptest_Eamp_stawise = nan(Nstas_dtstar,Ntest);
    a2ptest_Ephi_stawise = nan(Nstas_dtstar,Ntest);
    for is = 1:Nstas_dtstar
        a2ptest_Eamp_stawise(is,:) = nanmean(a2ptest_Eamp(any(is1is2==is,2),:));
        a2ptest_Ephi_stawise(is,:) = nanmean(a2ptest_Ephi(any(is1is2==is,2),:));
    end
    % baseline is each station's own middling amp and phi fits. 
    Eamp_baseline = (a2ptest_Eamp_stawise(:,test_a2pwt==parms.inv.amp2phiwt));
    Ephi_baseline = (a2ptest_Ephi_stawise(:,test_a2pwt==parms.inv.amp2phiwt));

    %% compare long-period amplitude to tstar. Shouldn't really be relationship
    

    %% some plots
    % if ifplot
        figure(31),clf;set(gcf,'pos',[400 213 703 973])
        for ixax = 1:5
            ax(ixax) = axes('Position',[0.1 1.01-(ixax*0.19) 0.85 0.16]);
            hold on;
        end
        for itest = 1:Ntest
            testcol = colour_get(log10(test_a2pwt(itest)),log10(max(test_a2pwt)),log10(min(test_a2pwt)),parula);
            err = a2ptest_Eamp_stawise(:,itest)+a2ptest_Ephi_stawise(:,itest);
            % tstar apparent 
            scatter(ax(1),1:Nstas_dtstar,a2ptest_dtstar(:,itest) - mean(a2ptest_dtstar,2),50./err,...
                'markerfacecolor',testcol,'linewidth',2,'MarkerEdgeColor','k')
            % dT apparent
            scatter(ax(2),1:Nstas_dtstar,a2ptest_dT(:,itest) - mean(a2ptest_dT,2),50./err,...
                'markerfacecolor',testcol,'linewidth',2,'MarkerEdgeColor','k')
            % A0 apparent
            scatter(ax(3),1:Nstas_dtstar,log(a2ptest_A0(:,itest)) - mean(log(a2ptest_A0),2),50./err,...
                'markerfacecolor',testcol,'linewidth',2,'MarkerEdgeColor','k')
            % mean station errors in Amp (^) and Phi (o) fits
    %         scatter([1:Nstas_dtstar]+(0.2*(itest/Ntest)-0.1),err,50,...
    %             'markerfacecolor',testcol,'linewidth',2,'MarkerEdgeColor','k')
            scatter(ax(4),[1:Nstas_dtstar]+(0.2*(itest/Ntest)-0.1),a2ptest_Eamp_stawise(:,itest)./Eamp_baseline,50,...
                'markerfacecolor',testcol,'linewidth',2,'MarkerEdgeColor','k','Marker','^')
            scatter(ax(4),[1:Nstas_dtstar]+(0.2*(itest/Ntest)-0.1),a2ptest_Ephi_stawise(:,itest)./Ephi_baseline,50,...
                'markerfacecolor',testcol,'linewidth',2,'MarkerEdgeColor','k','Marker','o')
        end 
    
        set(ax(5),'xtick',[1:Nstas_dtstar],'xticklabel',stas,'xlim',[0.8 (Nstas_dtstar+0.2)])
        set(ax(1:4),'xtick',[1:Nstas_dtstar],'xticklabel',[],'xlim',[0.8 (Nstas_dtstar+0.2)])
        ylabel(ax(1),'dtstar-variation','fontsize',16)
        ylabel(ax(2),'dT-variation','fontsize',16)
        ylabel(ax(3),'log(A0)-variation','fontsize',16)
        ylabel(ax(4),{'Error-variation'; '(relative to mean)'},'fontsize',16)
    % end

    %% station-wise metrics for possible focusing
 
    % absolute variation in peak-to-peak (amp-to-phi controlled) parameters
    maxrange_dtstar = abs(a2ptest_dtstar(:,iallamp) - a2ptest_dtstar(:,iallphi));
    maxrange_dT     = abs(a2ptest_dT(:,iallamp) - a2ptest_dT(:,iallphi));
    maxrange_dlgA0  = abs(log(a2ptest_A0(:,iallamp)) - log(a2ptest_A0(:,iallphi)));
    % worsening of fit to un-emphasized data types
    worsefit_amp = a2ptest_Eamp_stawise(:,iallphi)./Eamp_baseline;
    worsefit_phi = a2ptest_Ephi_stawise(:,iallamp)./Ephi_baseline;

    % focusing metric
    focusscore = prod([maxrange_dtstar,maxrange_dT,maxrange_dlgA0,worsefit_amp,worsefit_phi],2);
    focusscore(isnan(focusscore)) = inf;


    % add focusing metric to plot
    if ifplot
        plot(ax(5),[1:Nstas_dtstar]+(0.2*(itest/Ntest)-0.1),focusscore,'rx','linewidth',2)
        yline(focus_thresh,'--r')
        ylabel(ax(5),'focus score','fontsize',16)
    end


    % threshold focusscore = 
    gdofgd = focusscore < focus_thresh;

    %% some more metrics
    % correlation coefficient
    C2 = corrcoef(a2ptest_dtstar);
    if C2(iallphi,iallamp) < -0.5
        gdofgd = false(Nstas_dtstar,1);
    end

    figh = plot_focustest_pretty(stas,Amat,phimat,fmids,wtmat,test_a2pwt,a2ptest_dtstar,a2ptest_A0,a2ptest_Eamp_stawise,a2ptest_Ephi_stawise,focusscore,focus_thresh);







end
