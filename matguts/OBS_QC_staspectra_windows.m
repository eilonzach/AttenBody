function [ goodwins ] = QC_staspectra_windows( Zraw,H1raw,H2raw,Praw,tt,Nwin,ifplot )
% [ goodwins ] = QC_staspectra_windows( Zraw,H1raw,H2raw,Praw,tt,Nwin,ifplot )
%  Function to go through the different noise windows and do a QC on which
%  ones to keep based on if any of the windows particularly increase the
%  standard deviation of the average spectrum.
%  
% Z. Eilon 03/2016

pb = [0.002 2]; %pass-band, in Hz
tolerance = 1.5; % multiple of sigma to throw out if exceeded
a_val = 0.05; % for the f-test

%% prelims
tt = tt-min(tt);
dt = tt(2)-tt(1);
samprate = 1./dt;
npts = length(tt);
overlap = ceil((length(tt)*Nwin - length(Zraw))/(Nwin - 1));

%% compute spectra

[spect_Z,FZ] = spectrogram(Zraw,flat_hanning(tt,0.1*dt*npts),overlap,npts,samprate,'yaxis');
[spect_1] = spectrogram(H1raw,flat_hanning(tt,0.1*dt*npts),overlap,npts,samprate,'yaxis');
[spect_2] = spectrogram(H2raw,flat_hanning(tt,0.1*dt*npts),overlap,npts,samprate,'yaxis');
[spect_P] = spectrogram(Praw,flat_hanning(tt,0.1*dt*npts),overlap,npts,samprate,'yaxis');

%% weight so only frequencies within pb count 
% fwt
f_wt = ones(size(spect_Z,1),1);
f_wt(FZ<pb(1)) = 0;
f_wt(FZ>pb(2)) = 0;

%% get variability of spectra
a_Z = abs(spect_Z);
a_1 = abs(spect_1);
a_2 = abs(spect_2);
a_P = abs(spect_P);

la_Z = log10(a_Z).*(f_wt*ones(1,Nwin));
la_1 = log10(a_1).*(f_wt*ones(1,Nwin));
la_2 = log10(a_2).*(f_wt*ones(1,Nwin));
la_P = log10(a_P).*(f_wt*ones(1,Nwin));

sla_Z = zeros(size(a_Z));
sla_1 = zeros(size(a_1));
sla_2 = zeros(size(a_2));
sla_P = zeros(size(a_P));
for ii = 1:Nwin, 
    sla_Z(:,ii) = smooth(la_Z(:,ii),50); 
    sla_1(:,ii) = smooth(la_1(:,ii),50); 
    sla_2(:,ii) = smooth(la_2(:,ii),50); 
    sla_P(:,ii) = smooth(la_P(:,ii),50); 
end

dsla_Z = sla_Z-mean(sla_Z,2)*ones(1,Nwin);
dsla_1 = sla_1-mean(sla_1,2)*ones(1,Nwin);
dsla_2 = sla_2-mean(sla_2,2)*ones(1,Nwin);
dsla_P = sla_P-mean(sla_P,2)*ones(1,Nwin);

if ifplot 
figure(97),clf
subplot(411), semilogx(FZ,sla_Z)
subplot(412), semilogx(FZ,sla_1)
subplot(413), semilogx(FZ,sla_2)
subplot(414), semilogx(FZ,sla_P)

figure(98),clf
subplot(411), semilogx(FZ,dsla_Z)
subplot(412), semilogx(FZ,dsla_1)
subplot(413), semilogx(FZ,dsla_2)
subplot(414), semilogx(FZ,dsla_P)
end % ifplot

%% cycle through trying to kill high-std-norm windows
moveon = false;
goodwins = [1:Nwin]';
while moveon==false
    clear('normvarZ','normvar1','normvar2','normvarP');
    for ii = 1:length(goodwins)
        ind = goodwins; ind(ii) = [];

        normvarZ(ii,1) = norm(std(dsla_Z(:,ind),0,2));
        normvar1(ii,1) = norm(std(dsla_1(:,ind),0,2));
        normvar2(ii,1) = norm(std(dsla_2(:,ind),0,2));
        normvarP(ii,1) = norm(std(dsla_P(:,ind),0,2));

    end
    ubernorm = [normvarZ,normvar1,normvar2,normvarP];
    penalty = ones(length(goodwins),1)*median(ubernorm) - ubernorm; % large penalty if norm is v. different from the median 
    penalty = sum(penalty,2); % sum penalty across components
    
%     figure(2), clf
%     plot([1:Nwin],detrend(ubernorm,'constant'),'o-')
%     figure(3), clf
%     plot([1:Nwin],sum(ubernorm,2),'o-')

    kill = penalty>tolerance*std(penalty);
    if isempty(kill); break; end
    trypenalty = penalty(~kill,:);
    
    if ftest(penalty,1,trypenalty,1) < a_val
        goodwins = goodwins(~kill);
        Nwin = length(goodwins);
        moveon = false; % loop back
    else
        moveon=true; % no reason to remove more windows
    end
end

if ifplot
figure(99),clf
subplot(411), semilogx(FZ,dsla_Z(:,goodwins))
subplot(412), semilogx(FZ,dsla_1(:,goodwins))
subplot(413), semilogx(FZ,dsla_2(:,goodwins))
subplot(414), semilogx(FZ,dsla_P(:,goodwins))
end


end

