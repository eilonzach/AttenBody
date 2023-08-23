function [odat,Dwt] = attenuate( idat,samprate,dtstar,dT,A,alp )
% odat = attenuate( idat,samprate,dtstar,dT,A,alp=0 )
%   
% function to attenuate and time-shift one trace using the zero-phase
% differential attenuation operator.
% 
% INPUTS:
%   idat      - input data waveform, [npts x 1].
%   samprate  - sample rate for data
%   dtstar    - if positive it disperses and diminishes the input waveform
%   dT        - if psitive it delays the input waveform
%   A         - station gain term - multiple of original
%   alp       - frequency dependence of attenuation
% 
% OUTPUTS:
%   odat      - output data waveform, [npts x 1]
% 
% NB the input waveform should be tapered with a 10% Hanning taper. The
% output will be similarly tapered.

if nargin<6 || isempty(alp)
    alp=0;
end

%% pad data
idat = idat(:); % make column
npts = length(idat);
idat = [zeros(1000,1);idat;zeros(1000,1)]; % pad

%% Atten + dT operator (f-domain)
[IDAT,ff] = fft_ze(idat,1./samprate);
w = 2*pi*ff(:);
wabs = abs(w);
w0 = 2*pi;

if alp==0
    Aw = exp(-0.5*wabs*dtstar);
    phiw = (1/pi)*dtstar*log(2*pi*exp(pi)./wabs);
    phiw(1) = 2.5*phiw(2)-phiw(3);% fudge to get first element
elseif alp>0
    Aw = exp( -0.5*dtstar.*wabs*(wabs/w0).^(-alp) );
    phiw = 0.5 * (wabs/w0).^(-alp) * dtstar * cot(alp*pi/2);
    phiw(1) = 2.5*phiw(2)-phiw(3);% fudge to get first element
end
    
Dwt = Aw .* exp(-1i*w.*phiw(:)); % make f-domain attenuation operator
     
%% Attenuate
qdat = real(ifft(IDAT.*Dwt));

%% Delay
sampshft = round(dT*samprate);
maxpts = length(qdat);
if sampshft>0
qdat = [zeros(sampshft,1);qdat(1:maxpts-sampshft)];
elseif sampshft<0
qdat = [qdat(1-sampshft:maxpts);zeros(-sampshft,1)];
end

%% Gain
odat = A*qdat;

%% taper at the end
odat = odat(1000+[1:npts]); % un-pad
odat = flat_hanning_win([1:npts],odat,1,npts,round(npts/10)); % taper


end

