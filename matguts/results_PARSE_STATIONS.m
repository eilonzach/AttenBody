function [ stas_parsed,vals_parsed,slats_parsed,slons_parsed,selevs_parsed,indgd ] = results_PARSE_STATIONS( stas,vals,slats,slons,selevs )
%[ vals_parsed,slats_parsed,slons_parsed,selevs_parsed ] = results_PARSE_STATIONS( stas,vals,slats,slons,selevs )
%
% Function to parse data from all my stations - removing any stations for
% which there are no results/vals, and collating stations that are in the
% same locations because of shifting OBS deplouments (e.g. J35A and J35C)
%
% N.B. Measurements at the 'two' stations will be nan-meaned (i.e. asuming
% that they are NaNs if there was no measurement for that val.)

Nin = length(stas);

%% subset to only stations that have vals
indgd = [1:Nin]';
indgd(isnan(nanmean(vals,2))) = [];

stas_parsed = stas(indgd);
vals_parsed = vals(indgd,:);
slats_parsed = slats(indgd);
slons_parsed = slons(indgd);
selevs_parsed = selevs(indgd);
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
            % we have a match!
%             fprintf('%s %s\n',stas_parsed{is1},stas_parsed{is2})
            doubles = [doubles;is1,is2];
            done = [done;is2];
            end
        end
    end
end
Nrep = length(doubles);

%% update the primary stations' info
for is = 1:Nrep
    ind1 = doubles(is,1);
    ind2 = doubles(is,2);
    
    stas_parsed{ind1} = [stas_parsed{ind1},'_',stas_parsed{ind2}(end)];
    slons_parsed(ind1) = nanmean([slons_parsed(ind1),slons_parsed(ind2)]);
    slats_parsed(ind1) = nanmean([slats_parsed(ind1),slats_parsed(ind2)]);
    selevs_parsed(ind1) = nanmean([selevs_parsed(ind1),selevs_parsed(ind2)]);
    vals_parsed(ind1,:) = nanmean([vals_parsed(ind1,:);vals_parsed(ind2,:)]); % nanmean over the two rows of vals (if vals is a column vector, will just nanmean the two scalars making up the two rows
end


%% wipe the secondary stations' info
if ~isempty(doubles)
    kill = doubles(:,2);
    stas_parsed(kill) = [];
    slons_parsed(kill) = [];
    slats_parsed(kill) = [];
    selevs_parsed(kill) = [];
    vals_parsed(kill,:) = [];
    indgd(kill,:) = [];
end

end