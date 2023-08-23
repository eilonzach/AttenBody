function [ resp,faxis ] = mkresp_evalresp( respdir,net,sta,chan,T,N,ev_yyyy_jjj )
% [ resp,faxis ] = mkresp_evalresp( respdir,net,sta,chan,T,N,ev_yyyy_jjj )
%   get response function for a given network, station, channel, length of
%   data window (T), and number of samples (N)

wd = pwd;

% make frequency vector ends, nfreqs
fmin = (1/T);
if mod(N,2)
     fmax   = ((N-1)/2)*(1/T);
     nfreqs = (N-1)/2;
else
     fmax = (N/2)*(1/T);
     nfreqs = N/2;
end     
% go to respdir and use evalresp tool
cd(respdir)
system(sprintf('/Users/zeilon/Work/Codes/evalresp-3.3.3/evalresp %s %s %s %f %f %.0f -u dis -s lin ',...
    sta,chan,ev_yyyy_jjj,fmin,fmax,nfreqs));

% read newly made AMP and PHASE files
ampfile = dir(['AMP.',net,'.',sta,'.*.',chan ]);
phafile = dir(['PHASE.',net,'.',sta,'.*.',chan ]);
A = dlmread(ampfile(1).name);
P = dlmread(phafile(1).name);

% parse frequency axis, amplitude and phase of response
if mod(N,2)
     faxis  = [0;A(:,1);-flip(A(1:end  ,1))];
     resp_a = [0;A(:,2); flip(A(1:end  ,2))];
     resp_p = [0;P(:,2);-flip(P(1:end  ,2))];
else
     faxis  = [0;A(:,1);-flip(A(1:end-1,1))];
     resp_a = [0;A(:,2); flip(A(1:end-1,2))];
     resp_p = [0;P(:,2);-flip(P(1:end-1,2))];
end 
resp = resp_a.*exp(1i*d2r(resp_p));

cd(wd)
end

