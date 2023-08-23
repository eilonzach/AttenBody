function [ trace ] = fixtrace( trace )
% function to fix various predictable problems with the trace obtained
% using irisFetch, including repeated channels, split datastreams for same
% channel, empty channels, etc.

%% fix extra channels/locations/split channels
ntr = length(trace);
del = [];
chans = unique({trace.channel}); % find all channels
for ic = 1:length(chans)         % loop over channels
    indc = find(strcmp(chans{ic},{trace.channel})); 
    
    % check for simple duplicates
    if length(indc)>1
        while length(indc)>1 && isequal(trace(indc(1)),trace(indc(2))) % keep removing duplicates while they exist
            trace(indc(2))=[]; % remove duplicate
            indc = find(strcmp(chans{ic},{trace.channel})); 
        end
        if length(indc)==1, continue; end % ony one trace for this channel; was a duplication problem - move on. 
    end
            
    
    locs = unique({trace(indc).location}); % find all locations
    sr = zeros(length(locs),1);
    for ilo = 1:length(locs)               % loop over location
        indlo = find(strcmp(locs{ilo},{trace(indc).location}));
        sc = 0;
        dat= [];
        tt = [];
        for ilc = 1:length(indlo)        % loop over unique channel and location
            indclo = indc(indlo(ilc));   % index of this channal and location
            if isempty(trace(indclo).data), del = [del,indclo]; continue; end % del and continue if empty
            sc = sc + trace(indclo).sampleCount;
            dat= [dat; trace(indclo).data]; % concat data from split channels
            tt = [tt ; linspace(trace(indclo).startTime,trace(indclo).endTime,trace(indclo).sampleCount)'];
            if ilc > 1, del = [del,indclo]; end
        end
        trace(indc(indlo(1))).sampleCount = sc;
        trace(indc(indlo(1))).data = dat;
        trace(indc(indlo(1))).tt = tt; % add time field
        sr(ilo) = trace(indc(indlo(1))).sampleRate;
    end
    if length(locs)>1
        del = [del, indc(~strcmp(locs(sr==max(sr)),{trace(indc).location}))]; % add to del traces with loc not at max samprate    
    end
end
del = unique(del);
% delete extra channels
trace(del) = [];

% %% fix non-matching data lengths, padding with zeros
% chans = unique({trace.channel}); % find all channels
% 
% for ic = 1:length(chans)         % loop over channels
%     t0 = min


end

