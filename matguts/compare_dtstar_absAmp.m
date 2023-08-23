% script to compare the values of dtstar that we're getting out of the
% combspectra method to the absolute amplitudes of long-period waves in the
% dataset. If some are truly more attenuated (high dtstar) then they should
% have smaller absolute amplitudes, too.
% 
% Designed to be run within/after the calc_3_combspectra script

Tmin = 1;
Tmax = 3;
Nstas = size(all_dat0_gd,2);

tt = [datwind(1):1./resamprate:datwind(2)-1./resamprate]';
datf = filt_quick(all_dat0_gd,1/Tmax,1/Tmin,1./resamprate);

datwf = zeros(size(datf));
for is = 1:Nstas
datwf(:,is) = flat_hanning_win(tt,datf(:,is),-10,40,5);
end

amprms = rms(datwf);

figure(4), clf
subplot(211)
plot(tt,datf)
subplot(212)
hold on
for is = 1:Nstas
plot(tt,datwf(:,is),'color',colour_get(amprms(is),max(amprms),min(amprms)))
end

figure(5), clf, hold on
plot(delta_tstar_pref,amprms/mean(amprms),'o')
plot(delta_tstar_pref(find(isob(indgd))),amprms(find(isob(indgd)))/mean(amprms),'o','MarkerFaceColor','b')
plot(delta_tstar_pref,max(abs(datwf))/mean(abs(max(datwf))),'or')
plot(delta_tstar_pref(find(isob(indgd))),max(abs(datwf(:,find(isob(indgd)))))/mean(abs(max(datwf))),'or','MarkerFaceColor','r')
xlabel('$\Delta t^*$','interpreter','latex','FontSize',22)
ylabel('Amplitude','interpreter','latex','FontSize',22)
title(sprintf('Filter: %.1f to %.1f second',Tmin,Tmax),'interpreter','latex','FontSize',22)
set(gca,'FontSize',14,'box','on')



pause