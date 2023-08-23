function [ As,phis,wts ] = run_comb( dat1,dat2,fltinfo,wdo1,wdo2,jbds,dt,pretime,maxphi,cor_c_skip,fmin_cskip,ifplot )
% [ As,phis,wts ] = run_comb( dat1,dat2,fltdat,wdo1,wdo2,jbds,dt,pretime,maxphi,corcskip,fmin_cskip,ifplot )
%  function to run through the comb of filters and calcultate the phase
%  shift and amplitude scaling (as well as weight assigned) for each
%  narrow-band filter. Option whether or not to use low frequencies to do
%  the cycle-skip correction - at frequencies above fmin_cskip

    fltinfo = fltinfo(:);
    Nwds = length(fltinfo);
    if cor_c_skip
        if unique(diff([fltinfo.fmid]')<0)
%             fprintf('filter freqs are descending in value... flipping\n')
            fltinfo = flipud(fltinfo);
            flipped_fs = true;
        else
            flipped_fs = false;
        end
    else
        flipped_fs = false;
    end

    fmids = [fltinfo.fmid]';
    ttws = jbds.*dt - pretime;
    % taper
    datt1 = dat1.*wdo1; % taper whole data series
    datt2 = dat2.*wdo1; % taper whole data series

    As = zeros(Nwds,1);
    phis = zeros(Nwds,1);
    wts = zeros(Nwds,1);
    %% loop over frequency bins
    for iw = 1:Nwds
        % pad with plenty of zeros for the filter
        datp1 = [zeros(1000,1);datt1;zeros(1000,1)]; 
        datp2 = [zeros(1000,1);datt2;zeros(1000,1)]; 

        % FILTER zerophase bandpass filter & window
        datpf1=filtfilt(fltinfo(iw).bb, fltinfo(iw).aa, datp1);
        datpf2=filtfilt(fltinfo(iw).bb, fltinfo(iw).aa, datp2);
        % lop off padding
        datf1 = datpf1(1001:end-1000);
        datf2 = datpf2(1001:end-1000);
        % window
        datf1 = datf1.*wdo2;
        datf2 = datf2.*wdo2;
        qdatwf1 = datf1(jbds);
        qdatwf2 = datf2(jbds);    

        %% things to fit the phase bearing in mind lower frequencies...

        % find observed phase shift
        [dcor]=xcortimes([qdatwf1,qdatwf2], dt, pretime, maxphi,0);
        phi_f_obs = diff(dcor);

            if cor_c_skip
                if fmids(iw) > fmin_cskip
                    Npi = -2:1:2; % try cycle shifts of up to 4pi in either direction
                    E = zeros(length(Npi),1);
                    phshifts = Npi*(1./fmids(iw));
                    for ip = 1:length(Npi)
                        lnfs = log(fmids(1:iw));
                        tryphis = [phis(1:iw-1);phi_f_obs+phshifts(ip)];
                        % least squares fit - error only in y
                        G = [lnfs,ones(size(lnfs))];
                        mobs = [G'*G]\G'*tryphis;
                        m = mobs(1); b = mobs(2);
        %                 if ifplot
        %                     figure(1), clf, hold on
        %                     plot(lnfs,tryphis,'o')
        %                     plot(lnfs(iw),tryphis(iw),'or')
        %                     plot(lnfs,m*lnfs + b,'g')
        % %                     pause(0.01)
        %                 end
                        e = [m*lnfs+b]-tryphis;
                        E(ip) = norm(e,2);
                    end
                
                phi_f_obs_prelim = phshifts(E==min(E)) + phi_f_obs;
                % in unlikely event that error is same either way, use one
                % that is closest to previous time shift
                if length(phi_f_obs_prelim)>1
                    phi_f_obs_prelim = phi_f_obs_prelim(mindex(phi_f_obs_prelim,phis(iw-1)));
                end
                % in super unlikely event that also both equally close to
                % prev time shift... use one that is closest to zero
                if length(phi_f_obs_prelim)>1
                    phi_f_obs_prelim = phi_f_obs_prelim(mindex(phi_f_obs_prelim));
                end

                % shift once to correct to nearest cycle
                qdatwf2s = interp1(ttws-phi_f_obs_prelim,qdatwf2,ttws,'linear',0)';
                % now re-do xcortimes to get best fit within a half cycle of that..
                [dcor]=xcortimes([qdatwf1,qdatwf2s], dt, pretime, .5./fmids(iw),0);
                phi_f_obs = diff(dcor) + phi_f_obs_prelim;
                end
            end
        
        % make phase-corrected time series
        qdatwf2s = interp1(ttws-phi_f_obs,qdatwf2,ttws,'linear',0)';

        % calc. observed amplitude factor
        A_f_obs = (qdatwf2s'*qdatwf1)/(qdatwf1'*qdatwf1);
        %     A_f_obs = 1./((datwf'*qdatwfs)/(datwf'*datwf)); %<< could do it other way too

        % make amplitude-corrected time series
        qdatwf2sa = qdatwf2s./A_f_obs;
        
        
        %% A-OK plot and record result at this freq
        if ifplot
            figure(2), clf, set(gcf,'pos',[30 150 1200,400]), hold on
            plot(ttws,qdatwf1,'k','LineWidth',2)
            plot(ttws,qdatwf2,'b','LineWidth',1.5)
            plot(ttws,qdatwf2s,'r','LineWidth',1.5)
            plot(ttws,qdatwf2sa,'g','LineWidth',1)
            title(sprintf('Period = %.2f, A$_{obs}$ = %.2f, $\\phi_{obs}$ = %.2f',...
            1./fmids(iw),A_f_obs,phi_f_obs),'Fontsize',18,'interpreter','latex')
%             xlim([-20 20])
             pause 
        end
       
        As(iw) = A_f_obs;
        phis(iw) = phi_f_obs;

        acor = (qdatwf1'*qdatwf2sa).^2./(qdatwf1'*qdatwf1)./(qdatwf2sa'*qdatwf2sa);
        wts(iw) = acor.^2;
    end %  loop on narrow-band windows

    if flipped_fs % flip back...
        As = flipud(As);
        phis = flipud(phis);
        wts = flipud(wts);
    end
    
end

