function [ cmap_out ] = cmap_condense( cmap_in )
%[ cmap_out ] = cmap_condense( cmap_in )
%   function to condense the middle of a colourmap to emphasise small
%   variations

[M,N] = size(cmap_in);

x = linspace(-M/2,M/2,M);
xx = erf(4*x/M)*M/2;
% plot(x,xx,'r',x,x,'b')

cmap_out = interp1(x,cmap_in,xx);


end

