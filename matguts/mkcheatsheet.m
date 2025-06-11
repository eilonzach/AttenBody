function phasedb = mkcheatsheet(cheatsheetphases,gcarc_range,evdep)
% phasedb = mkcheatsheet(cheatsheetphases,gcarc_range,evdep)
%   Function to make a table of arrival times of relevant phases to plot on
%   top of the record section, to alert user to phase interference

% add depth phases, if relevant
if evdep
    cheatsheetphases = [cheatsheetphases,'pP','sS','sSKS','pPP'];
end
Nph = length(cheatsheetphases);
% now put into comma-separated string for taup
phs = cheatsheetphases{1};
for ip = 2:Nph, phs = [phs,',',cheatsheetphases{ip}]; end

% make array of distances
gc = [round(min(gcarc_range),1)-1 : 0.5 : round(max(gcarc_range),1)+1]'; 
Ngc = length(gc);

% prep output, and use tauptime
phasedb = struct('phases',{cheatsheetphases},'gcdists',gc,'phasetimes',nan(Ngc,Nph),'Nph',Nph,'Ngc',Ngc);
for id = 1:Ngc
    tpt = tauptime('p',phs,'z',evdep,'d',gc(id));
    % make phases unique
    [~,unqph] = unique({tpt.phase},'stable');
    tpt = tpt(unqph);
    % now associate found phases to desired phases
    [~,~,indP] = intersect({tpt.phase},cheatsheetphases,'stable');
    phasedb.phasetimes(id,indP) = [tpt.time];
end

%% While testing, plot
% figure(453),clf
% plot(phasedb.phasetimes,phasedb.gcdists,'--k')
% hold on
% for ip = 1:phasedb.Nph
%     gcplt = phasedb.gcdists(randi(phasedb.Ngc));
%     ttplt = interp1(phasedb.gcdists,phasedb.phasetimes(:,ip),gcplt);
%     text(ttplt+1,gcplt,phasedb.phases{ip},'color','r')
% end
% set(gca,'ydir','reverse')



end

