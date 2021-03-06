% function [paramf]=filtparam(param,lis)
% function qui filtre le param?tre param avec
% un filtre butterworth.
% 1/(lis/2) = fr?qence de coupure

function [paramf]=filtparam(param,lis)

[b,a]=butter(2,2/lis);
paramf=NaN*ones(size(param));
[nsta ndep]=size(param);
for ista=1:nsta
  ifin=isfinite(param(ista,:));
  inan=find(ifin==1);
  if length(inan)>lis/2
    paramf(ista,ifin)=filtfilt(b,a,param(ista,ifin));
  end
end

