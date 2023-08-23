function [ cmap_out ] = cmap_stretch( cmap_in )
%[ cmap_out ] = cmap_stretch( cmap_in )
%   function to stretch the middle of a colourmap to mute small variations

[M,N] = size(cmap_in);

x = linspace(-M/2,M/2,M);
xx = sign(x).*x.^2./(M/2);
% plot(x,xx,'r',x,x,'b')

cmap_out = interp1(x,cmap_in,xx);


end

