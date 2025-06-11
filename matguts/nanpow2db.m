function nan_dB = nanpow2db(y)
% Convert power to dB and turn bad values to nan

% We want to guarantee that the result is an integer
% if y is a negative power of 10.  To do so, we force
% some rounding of precision by adding 300-300.
nan_dB = (10.*log10(y)+300)-300;
nan_dB(y(:)<=0) = nan;
end