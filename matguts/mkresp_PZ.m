function [ resp,faxis ] = mkresp_PZ( pp,zz,gain,T,N )
% [ resp,faxis ] = mkresp_PZ( pp,zz,gain,T,N )
%   get response function for a given network, station, channel, length of
%   data window (T), and number of samples (N)


if mod(N,2)
     faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/T);
else
     faxis = [0:N/2,-N/2+1:-1]*(1/T);
end

w = faxis.*2*pi;

resp = ones(size(w));
for ip = 1:length(pp)
    resp = resp./(1i*w - pp(ip));
end
for ip = 1:length(zz)
    resp = resp.*(1i*w - zz(ip));
end

resp = resp*gain;