init_param

ncstartup

inpath='/home1/capchat/reference/';
fname='bathy_etopo2.nc';
f_cdf = netcdf([inpath fname],'read');


lat= f_cdf{'latitude'}(:);
lon= f_cdf{'longitude'}(:);

ilat=find(lat>=lat_min & lat<=lat_max);
ilon=find(lon>=lon_min & lon<=lon_max);