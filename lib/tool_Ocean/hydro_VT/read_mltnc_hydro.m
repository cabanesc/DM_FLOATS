function [lat,lon,juld,deph,pres,temp,tpot,psal,vorp,sig0]=read_mltnc_hydro(dirname,ficname)

% Initialisation
%inpath='/home/capchat/HYDRO_ATLN/MLT_NC/';

ncstartup

f_cdf = netcdf([dirname ficname],'read');

latb= f_cdf{'LATITUDE_BEGIN'}(:);
lonb= f_cdf{'LONGITUDE_BEGIN'}(:);
juldb=f_cdf{'JULD_BEGIN'}(:);
latb(latb==-9999)=NaN;
lonb(lonb==-9999)=NaN;
juldb(juldb==-9999)=NaN;

late= f_cdf{'LATITUDE_END'}(:);
lone= f_cdf{'LONGITUDE_END'}(:);
julde=f_cdf{'JULD_END'}(:);
late(late==-9999)=NaN;
lone(lone==-9999)=NaN;
julde(julde==-9999)=NaN;

lat=mynanmean([latb late]')';
lon=mynanmean([lonb lone]')';
juld=mynanmean([juldb julde]')';

clear latb lonb juldb late lone julde 

tabvar=['DEPH';
        'PRES';
	'TEMP';
	'TPOT';
	'PSAL';
	'VORP';
	'SIG0'];
           
for ivar=1:7
  tempo=f_cdf{tabvar(ivar,:)}(:);
  fillvalue=f_cdf{tabvar(ivar,:)}.FillValue_(:);
  tempo(find(tempo==fillvalue))=NaN;
  eval([lower(tabvar(ivar,:)) '=tempo;']);
end

close(f_cdf)