% function [paramf]=filt_param(param,lis)
% function qui filtre le paramètre param avec
% un filtre butterworth.
% 1/(lis/2) = fréqence de coupure

function [paramf]=filt_param(param,lis)

[b,a]=butter(2,2/lis);
paramf=NaN*ones(size(param));
[nsta ndep]=size(param);
for ista=1:nsta
  ifin=isfinite(param(ista,:));
  paramf(ista,ifin)=filtfilt(b,a,param(ista,ifin));
end

