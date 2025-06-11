%% Assign stations unique station No. - merging repeat/near-identical stas
clc

%% only real parameter - distance within which to exactly associate stations
drlim = 0.3./111; % = 300m; lateral distance to consider same station. In deg

% project details
dbname = 'EARdb';
dbdir = '/Users/zeilon/Dropbox/Work/EARdb/'; % include final slash

% dbname = 'FRES_PILOT';
% dbdir = '~/Dropbox/Work/FRES_PILOT/'; % include final slash


%% Preliminaries
wd = pwd;
addpath('matguts')
cd(dbdir);
run([dbdir,dbname,'_startup.m']);
% station details
load([infodir,'/stations'],'stainfo');

%% inter-sta distances
try 
    load([infodir,'/stations'],'stasta_dist_deg');
catch
    cd([dbdir,'INFO/'])
    run('mk_inter_sta_dist.m')
    try
        load([infodir,'/stations'],'stasta_dist_deg');
    cd(wd)
    catch
        cd(wd)
        error('no station-station distances')
    end
end


%% Assign unique station ID
% count_unqsta = 0;
staunqid = nan(stainfo.nstas,1);
for is = 1:stainfo.nstas
    % check if already assigned. If so, move on
    if ~isnan(staunqid(is))
        continue
    end
    fprintf('%2s-%s\n',stainfo.nwk{is},stainfo.stas{is})

    % new station, assign new id, corresponding to the row in the stainfo
    % matrix. So some ids (corresponding to repeat stations) will never be
    % assigned. The number of stations is not the highest unique id.
%     count_unqsta = count_unqsta+1;
    staunqid(is) = is;

    % find any obviously matching stations (same station and network)
    samenw = strcmp(stainfo.nwk,stainfo.nwk{is}); % find matching network
    samesta = regexpany(stainfo.stas,stainfo.stas{is}); % find matching station name string
    samenwsta = find(samenw & samesta);
    if length(samenwsta)>1
        staunqid(samenwsta(2:end)) = staunqid(is);
        [stainfo.nwk(samenwsta),stainfo.stas(samenwsta)]
    end

    % now find very nearby stations
    nearlalo = find(stasta_dist_deg(:,is) < drlim);
    nearlalo(nearlalo==is) = []; % not the same station
    nearlalo(~isnan(staunqid(nearlalo))) = []; % not already attributed
    if ~isempty(nearlalo)
%         staunqid(samenwsta(2:end)) = staunqid(is);
        fprintf('   %2s-%-5s %7.4f %7.4f\n',...
            stainfo.nwk{is},stainfo.stas{is},...
            stainfo.slats(is),stainfo.slons(is))
        for iiis = 1:length(nearlalo)
            fprintf('   %2s-%-5s %7.4f %7.4f  (%4.0f m)\n',...
                stainfo.nwk{nearlalo(iiis)},stainfo.stas{nearlalo(iiis)},...
                stainfo.slats(nearlalo(iiis)),stainfo.slons(nearlalo(iiis)),...
                stasta_dist_deg(is,nearlalo(iiis))*111e3)
        end


%         [stainfo.slats(is),stainfo.slons(is);stainfo.slats(nearlalo),stainfo.slons(nearlalo)]
%         [stainfo.nwk(nearlalo),stainfo.stas(nearlalo)]
        % call these the same station!
        staunqid(nearlalo) = staunqid(is);
    end 
end

stainfo.staunqid = staunqid;

%% save
save([infodir,'/stations'],'stainfo','-append');


  
