function [lat,lon,juld,deph,pres,param]=read_mltnc_hydro_1param(dirname,paramname)

% Initialisation
ncstartup

ldir=length(dirname);
inpath='/home1/penfret/HYDROCEAN/MLT_NC/';
if isempty(findstr(dirname,'ovid06')) == 0
    f_cdf = netcdf([dirname],'read');
else
    f_cdf = netcdf([inpath dirname],'read');
end

lat= f_cdf{'LATITUDE_BEGIN'}(:);
lon= f_cdf{'LONGITUDE_BEGIN'}(:);
juld=f_cdf{'JULD_BEGIN'}(:);

tabvar=['DEPH';
        'PRES';
	'TEMP';
	'TPOT';
	'PSAL';
	'VORP';
	'SIG0';
	'SIG1'];
      
varval=[1 2];
for ivar=3:8
  if strcmp(paramname,tabvar(ivar,:)) == 1
    varval=[varval ivar];
  end
end

for ivar=varval
  tempo=f_cdf{tabvar(ivar,:)}(:);
  fillvalue=f_cdf{tabvar(ivar,:)}.FillValue_(:);
  tempo(find(tempo==fillvalue))=NaN;
  if ivar == 1 | ivar == 2
    eval([lower(tabvar(ivar,:)) '=tempo;']);
  else
    param=tempo;
  end
end
