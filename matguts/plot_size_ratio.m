function [ height_to_width ] = plot_size_ratio( lonlims,latlims )
% [ height_to_width ] = plot_size_ratio( lonlims,latlims )
%   Function to calculate lat distance vs lon distance to have right
%   dimensions on plots

height = distance_km(latlims(1),lonlims(1),latlims(2),lonlims(1));
width  = distance_km(latlims(1),lonlims(1),latlims(1),lonlims(2));

height_to_width = height/width;


end

