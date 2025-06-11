function [ stas_parsed,vals_parsed,slats_parsed,slons_parsed,selevs_parsed,ondate_parsed,offdate_parsed,duration_parsed,indgd ] = results_PARSE_STATIONS( stas,vals,slats,slons,selevs,ondate,offdate )
%[ stas_parsed,vals_parsed,slats_parsed,slons_parsed,selevs_parsed,ondate_parsed,offdate_parsed,duration_parsed,indgd ] = results_PARSE_STATIONS( stas,vals,slats,slons,selevs,ondate,offdate )
%
% Function to parse data from all my stations - removing any stations for
% which there are no results/vals, and collating stations that are in the
% same locations because of shifting OBS deplouments (e.g. J35A and J35C)
%
% N.B. Measurements at the 'two' stations will be nan-meaned (i.e. asuming
% that they are NaNs if there was no measurement for that val.)

Nin = length(stas);

if nargin < 6
    ondate = zeros(Nin,1);
    offdate = zeros(Nin,1);
end



%% subset to only stations that have vals
indgd = [1:Nin]';
indgd(isnan(nanmean(vals,2))) = [];

stas_parsed = stas(indgd);
vals_parsed = vals(indgd,:);
slats_parsed = slats(indgd);
slons_parsed = slons(indgd);
selevs_parsed = selevs(indgd);
ondate_parsed = ondate(indgd);
offdate_parsed = offdate(indgd);
duration_parsed = offdate_parsed - ondate_parsed;
Ngd = length(indgd);

%% find stations that are coincident names + locations
doubles = [];
done = [];
for is1 = 1:Ngd
    if any(is1==done), continue; end
    for is2 = is1+1:Ngd
        if any(is2==done), continue; end
        namel = length(stas_parsed{is1}); if namel<4, continue; end % don't do for short sta names
        if strcmp(stas_parsed{is1}(1:end-1),stas_parsed{is2}(1:end-1)) %  1:N-1 characters in name match
            if distance(slats_parsed(is1),slons_parsed(is1),slats_parsed(is2),slons_parsed(is2)) < 0.1 % distance is less than 0.1 deg
            %we have a match!
            % fprintf('MATCH %s %s\n',stas_parsed{is1},stas_parsed{is2})
            doubles = [doubles;is1,is2];
            done = [done;is2];
            end
        end
    end
end
Nrep = size(doubles,1);

%% update the primary stations' info
for is = 1:Nrep
    ind1 = doubles(is,1);
    ind2 = doubles(is,2);
    
    stas_parsed{ind1} = [stas_parsed{ind1},'_',stas_parsed{ind2}(end)];
    slons_parsed(ind1) = nanmean([slons_parsed(ind1),slons_parsed(ind2)]);
    slats_parsed(ind1) = nanmean([slats_parsed(ind1),slats_parsed(ind2)]);
    selevs_parsed(ind1) = nanmean([selevs_parsed(ind1),selevs_parsed(ind2)]);
    vals_parsed(ind1,:) = nanmean([vals_parsed(ind1,:);vals_parsed(ind2,:)]); % nanmean over the two rows of vals (if vals is a column vector, will just nanmean the two scalars making up the two rows
    % calculate on duration - not so simple
    if ondate_parsed(ind1) < ondate_parsed(ind2) && offdate_parsed(ind1)  < ondate_parsed(ind2) % ind1 all before ind2
        duration_parsed(ind1) = duration_parsed(ind1)+duration_parsed(ind2);
    elseif ondate_parsed(ind2) < ondate_parsed(ind1) && offdate_parsed(ind2)  < ondate_parsed(ind1) % ind2 all before ind1 
        duration_parsed(ind1) = duration_parsed(ind1)+duration_parsed(ind2);
    else % overlapping 
        duration_parsed(ind1) = max([offdate_parsed(ind1),offdate_parsed(ind2)]) - min([ondate_parsed(ind1),ondate_parsed(ind2)]);
    end
    ondate_parsed(ind1) = min([ondate_parsed(ind1),ondate_parsed(ind2)]);
    offdate_parsed(ind1) = max([offdate_parsed(ind1),offdate_parsed(ind2)]);
end


%% wipe the secondary stations' info
if ~isempty(doubles)
    kill = doubles(:,2);
    stas_parsed(kill) = [];
    slons_parsed(kill) = [];
    slats_parsed(kill) = [];
    selevs_parsed(kill) = [];
    vals_parsed(kill,:) = [];
    ondate_parsed(kill) = [];
    offdate_parsed(kill) = [];
    duration_parsed(kill) = [];
    indgd(kill,:) = [];
end

end