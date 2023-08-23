function plot_filter_comb( flos,fhis,dt,npol,TorF)
% function plot_filter_comb( flos,fhis,dt,npol )
% 
% flos - vector of low frequency end of bandpass
% fhis - vector of high frequency end of bandpass
% dt - 1./samprate
% npol - number of poles for filter
if nargin<5
    TorF = 'T';
end

Tmids = 1./(0.5*(flos+fhis));
Nwds = length(flos);

fmin = 1e-3;
fmax = 1e2;

% make freq. axis
faxis = logspace(log10(fmin),log10(fmax),2000)';



figure(33), clf, set(gcf,'position',[360   514   900   800]);
hold on

for iw = 1:Nwds
% fprintf('Tlo = %.2f  Thi = %.2f\n',1./flos(iw),1./fhis(iw))
[zz,pp,gain]=butter(npol, [flos(iw), fhis(iw)].*dt*2);

SOS = zp2sos(zz,pp,gain);

[h,w] = freqz(SOS,faxis*dt*2*pi); 

h(1./faxis < 0.3) = 0;


if strcmp(TorF,'T')
    plot(1./(faxis),abs(h),'Linewidth',2,'color',colour_get(Tmids(iw),max(Tmids),min(Tmids))); 
elseif strcmp(TorF,'F')
    plot(faxis,abs(h),'Linewidth',2,'color',colour_get(Tmids(iw),max(Tmids),min(Tmids))); 
end

end % loop on windows

if strcmp(TorF,'T')
    xlim([0 30])
    xlabel('\textbf{Period (s)}','FontSize',24,'interpreter','latex')
elseif strcmp(TorF,'F')
    xlabel('\textbf{frequency (Hz)}','FontSize',24,'interpreter','latex')
    set(gca,'xscale','log','xlim',[1./40, 1.5])
end

ylabel('\textbf{Comb filter amplitude}','FontSize',24,'interpreter','latex')
grid on, box on
set(gca,'FontSize',18,'LineWidth',2)
ylim([0 1.1])

end

