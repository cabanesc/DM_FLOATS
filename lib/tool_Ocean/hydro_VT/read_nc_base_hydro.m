% 
% function [lat,lon,juld,platform,varargout]=read_nc_base_hydro(file_name,varargin)
%
% Programme de lecture des fichiers créer par le programme base_hydro.m de
% C. Lagadec
%
% Creation : V. Thierry, jan 2009
% 

function [lat,lon,juld,platform,varargout]=read_nc_base_hydro(file_name,varargin)

%ncstartup
f_cdf = netcdf(file_name,'read');

% Reads general informations

%data_type   = f_cdf{'DATA_TYPE'}(:);
%format_version = f_cdf{'FORMAT_VERSION'}(:);
platform  = f_cdf{'CRUISE_NAME'}(:);
%data_centre  = f_cdf{'DATA_CENTRE'}(:);
%data_creation  = f_cdf{'DATE_CREATION'}(:);
%data_update  = f_cdf{'DATE_UPDATE'}(:);
%cycnum=f_cdf{'CYCLE_NUMBER'}(:);

lat= f_cdf{'LATITUDE_BEGIN'}(:);
lon= f_cdf{'LONGITUDE_BEGIN'}(:);
juld=f_cdf{'JULD_BEGIN'}(:);

if nargin > 1
  for i=2:nargin
    tempo=varargin{i-1};
    varia=f_cdf{tempo}(:);
    fillvalue=f_cdf{tempo}.FillValue_(:);
    varia(find(varia==fillvalue))=NaN;
    varargout{i-1}=varia;
  end
end

close(f_cdf)
  


