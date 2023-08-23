function [ slab_contrs ] = McCrory_slab_depth
% function to pull out McCrory slab depth contours

ctrfile = '~/Work/CASCADIA/DATA/McCrory2012/dataset_S8_slabdepth_handcontour';

fid = fopen(ctrfile,'r');
A = textscan(fid,'%f %f %f','headerlines',1,'delimiter',',');
fclose(fid);

clevs = unique(A{3},'stable');

slab_contrs = struct([]);
for ic = 1:length(clevs)
    slab_contrs(ic).depth = clevs(ic);
    slab_contrs(ic).lat = A{2}(A{3}==clevs(ic));
    slab_contrs(ic).lon = A{1}(A{3}==clevs(ic));
end

