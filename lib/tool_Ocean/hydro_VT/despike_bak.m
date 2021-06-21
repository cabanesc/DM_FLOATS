
function dataout = despike_bak(data,s)

%
% function dataout = despike_bak(data,s)
%
% 5-point median despike
% scan the data with a 5-point window and
% set the value to NaN if it differs from the median value by more
% than the spike size, s.
% NaNs are discarded before applying the 5-point window.
%
% BAK 19 Nov 2008
%

dataout = nan+data;

n = length(data);
ki = 1:n;

% First remove any absent values from input
knan = find(isnan(data));
ki(knan) = [];
data(knan) = [];

% not central window for data at ends
k = 1;
while k <= length(data)
    if length(data) < 5
        fprintf(2,'%s\n','fewer than 5 good data remain')
        fprintf(2,'%s\n','terminating search for spikes')
        break
    end
    %     keep the good ones, throw out the spikes; make a note of which data cycles are kept
    if k < 3
        kw = 3;
    elseif  k > length(data)-2
        kw = length(data)-2;
    else
        kw = k;
    end
    %     keyboard
    d5 = data(kw-2:kw+2);
    s5 = sort(d5);
    if abs(d5(3+k-kw)-s5(3)) > s
        data(k) = [];
        ki(k) = [];
        continue
    else
        k = k+1;
        continue
    end
end

dataout(ki) = data;

return

