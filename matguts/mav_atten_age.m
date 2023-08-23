function [ dat_mav,age_mav ] = mav_atten_age( dat_obs,age_obs,awin,wt_obs )
% [ dat_mav,age_mav ] = mav_atten_age( dat_obs,age_obs,[awin=1] )

if nargin < 3 || isempty(awin)
    awin = 2;
end
if nargin < 4 || isempty(wt_obs)
    wt_obs = ones(size(dat_obs));
end
    
hwin = awin/2;

age_mav = [0:0.5:ceil(max(age_obs))]';
dat_mav = nan(size(age_mav));
for ia = 1:length(age_mav)
    mxa = age_mav(ia)+hwin;
    mna = age_mav(ia)-hwin;
    d = dat_obs(age_obs<=mxa & age_obs >=mna);
    w = wt_obs(age_obs<=mxa & age_obs >=mna);
    dat_mav(ia) = nanmean(d.*w)./nanmean(w);
end


end

