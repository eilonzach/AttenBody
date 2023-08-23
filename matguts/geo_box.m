function [ box ] = geo_box( w,h,a,c,opt )
%[ box ] = geo_box( w,h,a,c )
%   function to calculate the five (inc. repeat) corners of a box centred
%   at "c" (in [lat,lon]) and "w" degrees wide, "h" degrees high, and
%   rotated clockwise by an angle of "a" in the local coordinate system. 
% box has two 5-length columns of [lat,lon]
%  option 1 is to do the simple pythagorian rotation
%  option 2 is to use reckon etc.

if nargin<5
    opt = 1;
end

if opt==1
    boxx = [w/2,h/2; w/2,-h/2; -w/2,-h/2; -w/2,h/2; w/2,h/2];
    R = [cosd(a),sind(a);-sind(a),cosd(a)];
    boxxr = (R*boxx')';
    box = [boxxr(:,2),boxxr(:,1)] + ones(5,1)*[c(1),c(2)];
end


if opt==2
    [arclen,az] = distance(c(1),c(2),c(1)+h/2,c(2)+w/2);
    [trla,trlo] = reckon(c(1),c(2),arclen,az+a);
    [brla,brlo] = reckon(c(1),c(2),arclen,(180-az)+a);
    [blla,bllo] = reckon(c(1),c(2),arclen,180+az+a);
    [tlla,tllo] = reckon(c(1),c(2),arclen,a-az);
    box = [trla,trlo;brla,brlo;blla,bllo;tlla,tllo;trla,trlo];
end


end

