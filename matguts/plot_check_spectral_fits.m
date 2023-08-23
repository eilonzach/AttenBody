function plot_check_spectral_fits(Nstas)
% plot...
gdevts = unique(datinfo(:,3));
figure(44), clf, set(gcf,'pos',[ 440 -139 1430 797])
figure(45), clf, set(gcf,'pos',[ 440 -139 1430 797])
kk = 0;
while kk < 7*8
    % pick random pair of stations   
    rs = randsample(Nstas,2); stind1 = stinds(rs(1)); stind2 = stinds(rs(2));
%     % pick random event
%     re = randsample(length(gdevts),1); ie = gdevts(re);
    % find rows of data matrices where these are
    ii = find( min(datinfo(:,[1,2]),[],2)==min(stind1,stind2) & max(datinfo(:,[1,2]),[],2)==max(stind1,stind2) & datinfo(:,3)==ie);
    if isempty(ii), continue, end
    % grab differential values according to our results

    [Amat_pred,phimat_pred] = pred_Amat_phimat( dtstar_allin1(rs),dT(rs),10.^(lnA(rs)),fmids,alp );
    
    kk = kk + 1;
    
    figure(44), subplot(8,7,kk), hold on, 
    scatter(fmids,d_Amp(ii),50*Wts(ii) + 0.01,'or','filled'),
    plot(fmids,log(Amat_pred),'Linewidth',2), 
    title(sprintf('Stas %.0f - %.0f, ev%.0f',rs,re))

    figure(45), subplot(8,7,kk), hold on, 
    scatter(fmids,d_phi(ii),50*Wts(ii) + 0.01,'og','filled'),
    plot(fmids,phimat_pred,'Linewidth',2), 
    title(sprintf('Stas %.0f - %.0f, ev%.0f',rs,re))

end

end