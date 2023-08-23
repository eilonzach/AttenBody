% script to sequentially load the grids of pairwise spectra and solve each
% for dtstar assuming a given alpha (not necessarily zero
% 
% the format of the results file will be the same as for RESULTS_EXTRACT,
%  i.e. a big nstas x nevts matrix of values with nans where unavailable

clear all
close all

% project details
dbname = 'EARdb';
% dbdir = '/Volumes/LaCie/Granite_EastAfrica/EARdb/'; % include final slash
dbdir = '/Users/zeilon/Dropbox/Work/EARdb/'; % include final slash

spd = 24*60*60; % seconds per day - to convert time in seconds to serial time

ifsave = true;

Tmin = 1;
Tmax = 25;
Nwds = 20;

fmids = 1./logspace(log10(Tmin),log10(Tmax),Nwds)';
Nf = length(fmids);

phacomp = 'ST';

% alp = 0.27;
% alps = sort([0:0.1:0.6,0.27]');
alps = sort([0:0.1:0.6,0.27]');
alps = 0.27;

amp2phiwt = 2;

%% QC parms
fmax = 1/3;
fwt = double(fmids<=fmax);
maxdist = 6; % in degrees
minWtsum = 5;
ifonlyOBS = 0.5; % or do 0.5 to allow if one of two is obs




%% Preliminaries
wd = pwd;
addpath('matguts')
cd(dbdir);
run([dbdir,dbname,'_startup.m']);
pairdir = [resdir,'PAIRSPECS/'];


%% =================================================================== %%
%% ==========================  GET TO WORK  ========================== %%
%% =================================================================== %%


Ealp = zeros(size(alps));
for ia = 1:length(alps)
alp = alps(ia);

parms = struct('alp',alp,'amp2phiwt',amp2phiwt,'fmax',fmax,'maxdist_deg',maxdist,'minWtsum',minWtsum,'ifonlyOBS',ifonlyOBS);


%% starts
ws = 2*pi*fmids;
w0 = 2*pi;
if alp==0
    Ax = -0.5*ws; 
    Px = 1 - log(ws/2/pi)./pi;
else
    Ax = -0.5 * w0^alp * ws.^(1-alp);
    Px = 0.5 * cot(alp*pi/2) * (ws/w0).^(-alp);
end

% event details
load([infodir,'/events'],'evinfo'); 
% station details
load([infodir,'/stations'],'stainfo');

%% distances
% fprintf('Computing all inter-station distances\n');
% dd = zeros(stainfo.nstas);
% for is1 = 1:stainfo.nstas
% for is2 = is1+1:stainfo.nstas
%     dd(is1,is2) = distance(stainfo.slats(is1),stainfo.slons(is1),stainfo.slats(is2),stainfo.slons(is2));
% end
% end
% dd = dd+dd';
% save([infodir,'/stations'],'stainfo','dd');
load([infodir,'/stations'],'dd'); % distance in degrees


%% results structures
% max N of measurements (this is the number of pairs, not of amplitudes/phases measured
Npair_max = handshake(round(stainfo.nstas/10))*evinfo.norids; % (say on average 1/10 of stas are operational at any time - if needed, matrices will lengthen.

d_Amp = nan(Nf*Npair_max,1);
d_phi = nan(Nf*Npair_max,1);
datinfo = nan(Nf*Npair_max,3);
Wts  = nan(Nf*Npair_max,1); % inverse of the weights

% N will be Npair(=count)*Nf 
% M will be 3*stainfo.nstas for lnAmp, dtstar, dT (in that order)
% so size of G will be [Npair*Nf x 3*stainfo.nstas] 

G_Amp = spalloc(Nf*Npair_max,3*stainfo.nstas,4*Nf*Npair_max);
G_phi = spalloc(Nf*Npair_max,3*stainfo.nstas,4*Nf*Npair_max);

count = 0; % the count is the number of station pairs ==> length of G will be Nf*count
for ie = 1:evinfo.norids
    pairfile = sprintf('%s%.0f_pairspecs_%s.mat',pairdir,evinfo.orids(ie),phacomp);
    if exist(pairfile,'file')~=2, continue, end
    load(pairfile)
    fprintf('\nDoing orid %.0f\n',ie)
    
    if size(pairwise.As,2) ~= Nwds, fprintf('Stored results seem to be from run with different Nwds\n'), continue; end
	
    Amat = pairwise.As;
    phimat = pairwise.phis;
    wtmat = double(pairwise.inds).*pairwise.wts;
    
    evNstas = unhandshake(size(Amat,1));
    [~,evINDsta] = intersect(stainfo.stas,sts,'stable');

    evcount = 0;
    for is1 = 1:evNstas
    for is2 = is1+1:evNstas
        evcount = evcount+1;   % handshake # within this event
        
        sind1 = evINDsta(is1); % overall index of station 1
        sind2 = evINDsta(is2); % overall index of station 2
        
        % QCs
        if dd(sind1,sind2) > maxdist % don't include pair if too far apart
            continue
        end
        if sum(wtmat(evcount,:)' .* fwt) < minWtsum, fprintf('badwt '),continue, end

        count = count+1;       % overall count
        
        yind = [1:Nf] + (count-1)*Nf;
        
        % add elements of G matrices
        G_Amp(yind,sind1)         = -1;
        G_Amp(yind,sind2)         = 1;
        G_Amp(yind,stainfo.nstas+sind1)   = -Ax;
        G_Amp(yind,stainfo.nstas+sind2)   = Ax;

        G_phi(yind,stainfo.nstas+sind1)   = -Px;
        G_phi(yind,stainfo.nstas+sind2)   = Px;
        G_phi(yind,2*stainfo.nstas+sind1) = -1;
        G_phi(yind,2*stainfo.nstas+sind2) = 1;

        d_Amp(yind) = log(Amat(evcount,:))';
        d_phi(yind) = phimat(evcount,:)';
        Wts(yind)  = wtmat(evcount,:)' .* fwt;
        
        datinfo(yind,1) = sind1; datinfo(yind,2) = sind2; datinfo(yind,3) = ie;
    end % loop on sta2
    end % loop on sta1
end % loop on orids
%% delete nans
G_Amp = G_Amp(1:Nf*count,:);
G_phi = G_phi(1:Nf*count,:);
d_Amp = d_Amp(1:Nf*count,:);
d_phi = d_phi(1:Nf*count,:);
Wts  = Wts(1:Nf*count,:);
datinfo = datinfo(1:Nf*count,:);

%% find stations with data
fprintf('\nStripping out nans\n')
isdat_Amp = reshape(full(any(G_Amp))',stainfo.nstas,3)';
isdat_phi = reshape(full(any(G_phi))',stainfo.nstas,3)';
gdstas = any([isdat_Amp;isdat_phi])';


%% correct for repeat stations
stinds = find(gdstas);
Ngd = length(stinds);
doubles = [];
done = [];
for is1 = 1:Ngd
    sind1 = stinds(is1);
    if any(done==sind1), continue; end
    for is2 = is1+1:Ngd
        sind2 = stinds(is2);
        if any(done==sind2), continue; end
        namel = length(stainfo.stas{sind1}); if namel<4, continue; end % don't do for short sta names
        
        if strcmp(stainfo.stas{sind1}(1:end-1),stainfo.stas{sind2}(1:end-1)) %  1:N-1 characters in name match
            if dd(sind1,sind2) < 0.1 % distance is less than 0.1 deg
            % we have a match!
            doubles = [doubles;sind1,sind2];
            done = [done;sind2];
            end
        end
        
    end
end
% correct doubles stas in Gs
for is = 1:length(doubles)
    xind1 = doubles(is,1) + [0,stainfo.nstas,2*stainfo.nstas];
    xind2 = doubles(is,2) + [0,stainfo.nstas,2*stainfo.nstas];
    G_Amp(:,xind1) = G_Amp(:,xind1) + G_Amp(:,xind2);
    G_phi(:,xind1) = G_phi(:,xind1) + G_phi(:,xind2);
end
% wipe doubled stas in gdstas
gdstas(done) = false;
stinds = find(gdstas);

%% subset Gs to stations with gdstas
Ngd = sum(gdstas);
xind = find(gdstas)*[1 1 1] + ones(Ngd,1)*[0,stainfo.nstas,2*stainfo.nstas];
xind = xind(:)'; % these are the indices of all columns in Gs to keep

G_Amp_do = G_Amp(:,xind);
G_phi_do = G_phi(:,xind);

Nobs = full(sum(G_Amp_do(:,1:Ngd)~=0))'/Nf; % count N of obs per station, summing non-zeros down rows and dividing by Nf fmids
Nwobs = full(sum(sparse(1:Nf*count,1:Nf*count,Wts,Nf*count,Nf*count,Nf*count)*G_Amp_do(:,1:Ngd)~=0))'/Nf; % count weighted obs per station, multiplying by weights (some of which are zero) and then summing non-zeros down rows and dividing by Nf fmids

%% make final matrices
fprintf('Making matrices\n')
constraint_lnAmp  = sparse(1,[1:Ngd]      ,1,1,3*Ngd);
constraint_dtstar = sparse(1,[1:Ngd]+Ngd  ,1,1,3*Ngd);
constraint_dT     = sparse(1,[1:Ngd]+2*Ngd,1,1,3*Ngd);

d_all = [d_Amp;d_phi;0;0;0];
G_all = [G_Amp_do;G_phi_do;constraint_dtstar;constraint_dT;constraint_lnAmp];

w_all = [amp2phiwt*Wts;Wts;1;1;1];
N = length(w_all);
spdw_all = sparse(1:N,1:N,w_all,N,N,N);


dw = spdw_all*d_all;
Gw = spdw_all*G_all;

    
%% solve the least squares problem
fprintf('Solving weighted least squares\n')
m = lsqr(Gw,dw,1e-6,1000);

% parse results
lnA = m(1:Ngd);
dtstar_allin1 = m([1:Ngd] + Ngd);
dT = m([1:Ngd] + 2*Ngd); 

%% analyse error

ew = (Gw*m - dw);
Ealp(ia) = ew'*ew

end % loop on alphas

if length(alps)>1
    figure(27), clf
    plot(alps,Ealp,'o-','linewidth',2)
    set(gca,'fontsize',15)
    xlabel('alpha','fontsize',17);
    ylabel('Sum of errors','fontsize',17)
    save2pdf(27,sprintf('allin1_stav_dtstar_%s_E_vs_alp',phacomp),resdir);

    fprintf('===================\nAssume you want to stop here, as testing different alphas. No plotting/saving.\n')
    return
else

%% simple plot


run([plotdir,'map_parameters'])
figure(1); clf, set(gcf,'pos',[200 200 600 800])
subplot(4,1,1:3), hold on
scatter(stainfo.slons(stinds),stainfo.slats(stinds),Nwobs+0.001,dtstar_allin1,'filled')
colormap(parula), caxis([-1 1])
% [chgrdX,chgrdY] = meshgrid(linspace(lonlims(1),-123,50),linspace(40,latlims(2),60)); 
% [ chgrdage,chrons,F ] = jdf_crust_age(chgrdY,chgrdX);
% % chtick = unique(chrons.age(~isnan(chrons.age))); % plot chrons
% chtick = 1:12; % plot Ma
% for ic = 1:length(chtick)
%     ch = chtick(ic);
%     contour(chgrdX,chgrdY,chgrdage,[ch ch],'LineWidth',2,'linestyle','--',...
%         'color',colour_get(ch,max(chtick),min(chtick),hot))
% end
xlim([lonlims]), ylim(latlims)

% dist2rft
stainfo.sta_dist2rft(stinds,1) = dist2line( [MER.lon(1),MER.lat(1)],[MER.lon(2),MER.lat(2)],[stainfo.slons(stinds),stainfo.slats(stinds)])';

% caxis([-0.8 0.5])
subplot(4,1,4), cla, hold on
ind = Nwobs>40;
scatter(stainfo.sta_dist2rft(stinds(ind)),dtstar_allin1(ind),Nwobs(ind)/2+0.001,'r','filled')
grid on
% if any(doubles)
%     [~,ind,~] = intersect(stinds,doubles(:,1),'stable');
%     scatter(sages(stinds(ind)),dtstar_allin1(ind),Nwobs(ind)/2+0.001,'b')  
% end
xlim([-6 6]), ylim([-1.5 1.5])

%% collate results
results = struct('dtstar',dtstar_allin1,'dT',dT,'lnA',lnA,'Nobs',Nobs,'Nwobs',Nwobs,...
                 'stas',{stainfo.stas(stinds)},'slats',stainfo.slats(stinds),'slons',stainfo.slons(stinds),...
                 'selevs',stainfo.selevs(stinds),'sdist2rft',stainfo.sta_dist2rft(stinds),...
                 'parms',parms);

if ifsave
    resfile = sprintf('allin1_stav_dtstar_%s_alp%03.0f',phacomp,100*alp);
    fprintf('Saving!\n')
    save([resdir,resfile],'results')
end



end
return


%% check some fits

% plot...
gdevts = unique(datinfo(:,3));
figure(44), clf, set(gcf,'pos',[ 440 -139 1430 797])
figure(45), clf, set(gcf,'pos',[ 440 -139 1430 797])
kk = 0;
while kk < 7*8
    % pick random pair of stations   
    rs = randsample(Ngd,2); stind1 = stinds(rs(1)); stind2 = stinds(rs(2));
    % pick random event
    re = randsample(length(gdevts),1); ie = gdevts(re);
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



