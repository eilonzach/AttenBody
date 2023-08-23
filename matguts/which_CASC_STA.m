function [ sta ] = which_CASC_STA( sta_name,slat,slon )
%[ sta ] = which_CASC_STA( sta_name,slat,slon )
%   Script to rename stations to remove ambiguity about which bloody
%   station it is - for double-named stations like M08C, J06A...
% appends _O or _L for ocean or land, respectively
% 
% Written, resentfully, by Z. Eilon 04/2016

switch sta_name
    case 'G03A'
        if abs(slat - 40.0591) < 0.02 && abs(slon - -126.1625) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 45.3153) < 0.02 && abs(slon - -123.2811) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'G03D'
        if abs(slat - 40.0587) < 0.02 && abs(slon - -126.1612) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 45.2115) < 0.02 && abs(slon - -123.2641) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'G05D'
        if abs(slat - 40.1430) < 0.02 && abs(slon - -127.8257) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 45.2422) < 0.02 && abs(slon - -121.3167) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'J06A'
        if abs(slat - 43.2515) < 0.02 && abs(slon - -128.8010) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 43.2515) < 0.02 && abs(slon - -120.1528) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M01C'
        if abs(slat - 49.1504) < 0.02 && abs(slon - -126.7222) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.8473) < 0.02 && abs(slon - -124.1221) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M02C'
        if abs(slat - 48.3069) < 0.02 && abs(slon - -125.6012) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.3920) < 0.02 && abs(slon - -122.8538) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M03C'
        if abs(slat - 47.8884) < 0.02 && abs(slon - -126.1046) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.2742) < 0.02 && abs(slon - -122.1220) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M04C'
        if abs(slat - 45.5584) < 0.02 && abs(slon - -125.1923) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.7826) < 0.02 && abs(slon - -121.8393) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M05C'
        if abs(slat - 46.1735) < 0.02 && abs(slon - -124.9345) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.3593) < 0.02 && abs(slon - -121.1457) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M06C'
        if abs(slat - 45.5287) < 0.02 && abs(slon - -124.9261) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.2047) < 0.02 && abs(slon - -120.4772) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M07A'
        if abs(slat - 44.8987) < 0.02 && abs(slon - -125.1168) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.3884) < 0.02 && abs(slon - -119.1711) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M08A'
        if abs(slat - 44.1187) < 0.02 && abs(slon - -124.8953) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.4483) < 0.02 && abs(slon - -118.3792) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    
    %% pre-named L
    case 'G03AL'
        if abs(slat - 40.0591) < 0.02 && abs(slon - -126.1625) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 45.3153) < 0.02 && abs(slon - -123.2811) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'G03DL'
        if abs(slat - 40.0587) < 0.02 && abs(slon - -126.1612) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 45.2115) < 0.02 && abs(slon - -123.2641) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'G05DL'
        if abs(slat - 40.1430) < 0.02 && abs(slon - -127.8257) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 45.2422) < 0.02 && abs(slon - -121.3167) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'J06AL'
        if abs(slat - 43.2515) < 0.02 && abs(slon - -128.8010) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 43.2515) < 0.02 && abs(slon - -120.1528) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M01CL'
        if abs(slat - 49.1504) < 0.02 && abs(slon - -126.7222) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.8473) < 0.02 && abs(slon - -124.1221) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M02CL'
        if abs(slat - 48.3069) < 0.02 && abs(slon - -125.6012) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.3920) < 0.02 && abs(slon - -122.8538) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M03CL'
        if abs(slat - 47.8884) < 0.02 && abs(slon - -126.1046) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.2742) < 0.02 && abs(slon - -122.1220) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M04CL'
        if abs(slat - 45.5584) < 0.02 && abs(slon - -125.1923) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.7826) < 0.02 && abs(slon - -121.8393) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M05CL'
        if abs(slat - 46.1735) < 0.02 && abs(slon - -124.9345) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.3593) < 0.02 && abs(slon - -121.1457) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M06CL'
        if abs(slat - 45.5287) < 0.02 && abs(slon - -124.9261) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.2047) < 0.02 && abs(slon - -120.4772) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M07AL'
        if abs(slat - 44.8987) < 0.02 && abs(slon - -125.1168) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.3884) < 0.02 && abs(slon - -119.1711) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M08AL'
        if abs(slat - 44.1187) < 0.02 && abs(slon - -124.8953) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.4483) < 0.02 && abs(slon - -118.3792) < 0.02
            sta = [sta_name(1:4),'_L'];
        end

%% pre-named O
    case 'G03AO'
        if abs(slat - 40.0591) < 0.02 && abs(slon - -126.1625) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 45.3153) < 0.02 && abs(slon - -123.2811) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'G03DO'
        if abs(slat - 40.0587) < 0.02 && abs(slon - -126.1612) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 45.2115) < 0.02 && abs(slon - -123.2641) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'G05DO'
        if abs(slat - 40.1430) < 0.02 && abs(slon - -127.8257) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 45.2422) < 0.02 && abs(slon - -121.3167) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'J06AO'
        if abs(slat - 43.2515) < 0.02 && abs(slon - -128.8010) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 43.2515) < 0.02 && abs(slon - -120.1528) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M01CO'
        if abs(slat - 49.1504) < 0.02 && abs(slon - -126.7222) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.8473) < 0.02 && abs(slon - -124.1221) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M02CO'
        if abs(slat - 48.3069) < 0.02 && abs(slon - -125.6012) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.3920) < 0.02 && abs(slon - -122.8538) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M03CO'
        if abs(slat - 47.8884) < 0.02 && abs(slon - -126.1046) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.2742) < 0.02 && abs(slon - -122.1220) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M04CO'
        if abs(slat - 45.5584) < 0.02 && abs(slon - -125.1923) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.7826) < 0.02 && abs(slon - -121.8393) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M05CO'
        if abs(slat - 46.1735) < 0.02 && abs(slon - -124.9345) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.3593) < 0.02 && abs(slon - -121.1457) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M06CO'
        if abs(slat - 45.5287) < 0.02 && abs(slon - -124.9261) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.2047) < 0.02 && abs(slon - -120.4772) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M07AO'
        if abs(slat - 44.8987) < 0.02 && abs(slon - -125.1168) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.3884) < 0.02 && abs(slon - -119.1711) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
    case 'M08AO'
        if abs(slat - 44.1187) < 0.02 && abs(slon - -124.8953) < 0.02
            sta = [sta_name(1:4),'_O'];
        elseif abs(slat - 41.4483) < 0.02 && abs(slon - -118.3792) < 0.02
            sta = [sta_name(1:4),'_L'];
        end
        
    otherwise
        sta = sta_name;
end

end

