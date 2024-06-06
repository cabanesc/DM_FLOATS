function percentiles = myquantile(data,samplePoints)

data(isnan(data))=[];
if isempty(data)==0
percentiles = interp1(linspace(0, 1, length(data)), sort(data), samplePoints)
else
percentiles = 999;
end