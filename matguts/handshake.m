function [ handshakea ] = handshake( a )
%[ handshakea ] = handshake( a )
% calculate handshakes between "a" persons, i.e. a(a-1)/2
handshakea = (a.^2 - a)./2;

end

