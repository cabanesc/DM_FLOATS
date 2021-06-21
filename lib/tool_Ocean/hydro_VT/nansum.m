
function y = nansum(x,dim)

% version of nansum to avoid stats toolbox
%
% columns with no finite values returns nan.
% The matlab version returns zero if there are zero non-nan values
%
% this version looks for finite values, which are non-nan and non-inf
% arithmetic with inf may be unpredictable, eg
% inf + 1 = inf;
% inf - inf = nan;
%
% Brian King


if nargin==1,
     dim = min(find(size(x) > 1)); % default is first non-unity dimension
     if isempty(dim), dim = 1; end
end

ok = isfinite(x);
x(~ok) = 0;
y = sum(x,dim);

n = sum(ok,dim);
n(n==0) = NaN; % set nans if there are zero non-nan values included in the sum
n = 0*n; % elements are now zeros with nans; this is a ( 0 or nan) mask for columns with zero non-nan values
y = y + n; % mask out the columns that contain zero non-nan values

return

