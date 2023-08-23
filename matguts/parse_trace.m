function [ trace ] = parse_trace( trace )
%[ trace ] = parse_trace( trace )
%  Parse through and potentially clean the trace object spat out by the
%  irisFETCH tool

k = 1;
allchans = {trace.channel};
if length(unique(allchans))<length(allchans) % must have repeated chans - maybe need to stitch
    donechans = zeros(size(allchans));
    while any(donechans==0)
        dochan = find(donechans==0,1,'first');
        newtrace(k) = trace(dochan);
        samechan = find(strcmp(allchans{dochan},allchans));
        if length(samechan) > 1
            newtrace(k).endTime = trace(samechan(end)).endTime;
            newtrace(k).sampleCount = sum([trace(samechan).sampleCount]);
            datall = [];
            for ic = 1:length(samechan)
                datall = [datall; trace(samechan(ic)).data];
            end
            newtrace(k).data = datall;
        end  
        donechans(samechan) = 1;
        k = k+1;
    end
trace = newtrace;
end

% % end

