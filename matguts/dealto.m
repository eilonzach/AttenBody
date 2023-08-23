function outstruct = dealto(instruct,fieldname,vals_array)
%  outstruct = dealto(instruct,fieldname,vals_array)
%  function to deal elements of some array - be that a cell array or a
%  vector into a particular field of a structure
%  e.g. let's say I have a structure "eq" with fields .sta, .artime and I
%  have a cell array of station names I'd like to deal into the first
%  field, and a vector array of times I'd like to deal into the second
%  field.
% 
%  I could do this with
%   tmp = num2struct(artimes)
%   [eq(ind).artime] = deal(tmp{:})
%  this function basically obviates having the annoying intermediate step
%  in your code.
% 
%  INPUTS:
%   instruct   - structure into which you want to put vals
%   fieldname  - name of field within structure to put vals into (string)
%   vals_array - vector (cell or other) array of vals to put into structure
% 
%  N.B. instruct and vals_array must be the same size.

if length(vals_array)~=1;
    if length(instruct)~=length(vals_array)
        error('Structure and array of values must be same size\n')
    end
end

if ~iscell(vals_array)
    vals_array = num2cell(vals_array);
end

[instruct.(fieldname)] = deal(vals_array{:});
outstruct = instruct;

end