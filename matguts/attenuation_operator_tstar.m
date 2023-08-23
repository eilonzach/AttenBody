function qdat = attenuation_operator_tstar(dat,tstar,dt,alpha,A0,dT0,f0,fflag) % do the forward model (attenuate, scale, shift)
% qdat = attenuation_operator_tstar(dat,tstar,dt,alpha,A0,dT0,f0,fflag) % do the forward model (attenuate, scale, shift)
% 
% dat must be the data, in columns (i.e. a column vector for one time
% series, or many columns for multiple time series)
% 
% dt is the sample interval. dT0 is the time shift
% 
% positive A0 numbers indicate more gain (higher amplitude)
% positive dT0 numbers indicate later arrival (slower)
%
% fflag = 1 means "dat" is actually the [DAT,ff] (f-domain) 
% fflag = 0 means "dat" is waveforms (t-domain) and must be fft-ed

    if nargin < 8 || isempty(fflag)
        if any(imag(dat)~=0,'all')
            fflag = 1;
        else
            fflag = 0;
        end
    end
    if nargin < 7 || isempty(f0)
        f0 = 1; % default freq is 1 H
    end
    if nargin < 6 || isempty(dT0)
        dT0 = 0; % default no time shift also
    end
    if nargin < 5 || isempty(A0)
        A0 = 1; % default no gain difference also
    end
    if nargin < 4 || isempty(alpha)
        alpha = 0; % default no f dependence
    end    

    % get sizing right, in case multiple traces or atten values
    if any(size(dat)==1), dat = dat(:); end % columnise if necessary
    Nsamp = size(dat,1);
    Nat = length(tstar); % number of attenuation values
    Ntr = size(dat,2); % number of traces

    % take fft if needed
    if ~fflag
        [DAT,ff] = fft_ze(dat,dt);
    elseif fflag
        DAT = dat(:,[1:end-1]);
        ff  = dat(:,end);
        Ntr = Ntr-1;
    end
    
    % check sizing
    if Ntr > 1 && Nat>1 && Ntr~=Nat, error('Must have same number of traces as attenuation values, or one/other should be singleton'); end
    if Nat > 1 && length(A0)~=Nat, A0 = A0*ones(Nat,1); end % extend A0 across dtstar span 
    if Nat > 1 && length(dT0)~=Nat, dT0 = dT0*ones(Nat,1); end % extend dT0 across dtstar span
    if Nat > 1 % turn atten etc. into rows
        tstar = tstar(:)';
        A0 = A0(:)';
        dT0 = dT0(:)';
    end

    % angular freq
    w = 2*pi*ff;
    w0 = 2*pi*f0;
    winf = 2*pi*exp(pi); % note this assumes f0 is 1Hz
%     wgrid = w*ones(1,Nat); % repeat down rows, in case multiple attens;
    Pw = zeros(Nsamp,Nat);
    if alpha == 0 % frequency independent
        Aw = exp(-0.5*abs(w)*tstar);
        phw = (log(winf)-log(abs(w)))*(tstar/pi) + ones(Nsamp,1)*dT0;
%         phw(abs(w)>w0) = 0;
%         Pw = exp(-1i.*wgrid.*phw);
%         Pw = exp(-1i.*w*tstar./pi  - 1i.*w.*dT0);
        
    elseif alpha > 0
        gamma = 0.5*cot(alpha*pi/2);
        Aw = exp(-0.5*(w0.^alpha).*tstar.*abs(w).^(1-alpha));
        phw = tstar.*gamma.*(w0./abs(w)).^alpha + dT0 ;
    end

    for ii = 1:Nat
        Pw(:,ii) = exp(-1i.*w.*phw(:,ii));
        Pw(w==0,ii) = 1;
    end

%     Pw = exp(-1i.*wgrid.*phw);
%     if any(isnan(Pw))
%         if all(find(isnan(Pw)) == find(wgrid==0))
%               Pw(wgrid==0) = 1; 
%         else
%             error('some nans in Pw where f is not 0')
%         end 
%     end
    Dw = Aw.*Pw;
    
    % apply attenuation operator and ifft
    qdat = real(ifft(DAT.*Dw));

    % amplify (faster to loop here, actually)
%     qdat = (ones(Nsamp,1)*A0).*qdat;
    for itr = 1:Ntr
        qdat(:,itr) = A0(itr)*qdat(:,itr);
    end
end