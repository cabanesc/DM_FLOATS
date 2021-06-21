%-------------------------------------------------------------------------------
%
  function [msg_error] = ...
           NCW_OA_field(OA_filename, longitude, latitude, depth, jul_rel, ...
                         CONFIG, PARAM, FIELD)
%
%-------------------------------------------------------------------------------
%  
%   This function creates and write a NetCDF files for a gridded 3D parameter
%
%  OA_filename : Full name of file to be created
%  longitude:	 Vector of longitudes (its length defines the longitude 
%                dimension)
%  latitude:	 Vector of latitudes (its length defines the latitude 
%                dimension)
%  depth:	     Vector of depth (its length defines the depth dimension)
%  CONFIG:	     structure describing the configuration
%  PARAM: 	     parameter ('TEMP', 'PSAL', 'TPOT', 'SIG0')
%  FIELD:	     3D griddes field FIELD(depth, latitude, longitude) - 
%         	     can be temperature (PARAM = 'TEMP') or salinity 
%                (PARAM = 'PSAL')
%
%-------------------------------------------------------------------------------
%	
%  version:
%  --------
%  1.01   30/04/2003   E. Autret
%       Creation
%
%  1.02   19/04/2005  F. Gaillard
%       Update (comments and date réference)
%
%  1.03   01/07/2005  F. Gaillard
%       new parameters added (TPOT, SIG0), 
%       factorisation of variables definition
%       variable saved as float
%
%  1.04   19/01/2006  F. Gaillard
%
%  1.05   18/01/2007  F. Gaillard
%        add sigma as vertical coordinate
%       
%-------------------------------------------------------------------------------


%  Messages:
msg_file = ['OA_ncwrite_fl : Probleme d''ouverture du fichier :' OA_filename];

%  ---------------------------------------------------------------

% ------------------------
%    Creates the new file
% ------------------------
nc = netcdf(OA_filename, 'clobber');
if isempty(nc)
   msg_error = msg_file;
   return
else
   msg_error = 'ok';  
end


% ----------------------------
%    Defines the dimensions
% ----------------------------
%  Fixed dimensions
nc('time')   = length(jul_rel); 
if PARAM(1)=='Z'
    ZCOO = PARAM(2:end);
    PARSAV = 'DEPH';
else
    ZCOO = 'depth';    
    PARSAV = PARAM;
end
nc(ZCOO)        = length(depth);
nc('latitude')  = length(latitude);
nc('longitude') = length(longitude);


% ----------------------------
%    Defines global attributes
% ----------------------------

nc.CONVENTIONS     = ncchar('COARDS');
nc.producer_agengy = ncchar('IFREMER');
nc.project_name    = ncchar(CONFIG.PROJECT_NAME);
nc.creation_time   = ncchar(datestr(now,30));
nc.software_version= ncchar(CONFIG.SOFTWARE_VERSION);
nc.product_version = ncchar(CONFIG.PRODUCT_VERSION);
nc.data_set        = ncchar(CONFIG.DATA_SET);
nc.data_manager    = ncchar(CONFIG.DATA_MANAGER);  
nc.estimate_date  = ncchar(CONFIG.DATE_EST);
nc.south_latitude = ncchar(CONFIG.SOUTH_LAT);
nc.north_latitude = ncchar(CONFIG.NORTH_LAT);
nc.west_longitude = ncchar(CONFIG.WEST_LONG);
nc.east_longitude = ncchar(CONFIG.EAST_LONG);

% ----------------------------------------
%   Variables and attributes:
% ----------------------------------------
%
nc{'time'} = ncfloat('time');
nc{'time'}.units = ncchar(['days since 1950/01/01 UTC 00:00:00']);

nc{'latitude'} = ncfloat('latitude');
nc{'latitude'}.units = ncchar('degree_north');
nc{'latitude'}.valid_min = ncfloat(-90);
nc{'latitude'}.valid_max = ncfloat(90);

nc{'longitude'} = ncfloat('longitude');
nc{'longitude'}.units = ncchar('degree_east');
nc{'longitude'}.valid_min = ncfloat(-180);
nc{'longitude'}.valid_max = ncfloat(180);


if strcmp(ZCOO,'depth')       
    nc{ZCOO} = ncshort('depth');
    nc{ZCOO}.units =ncchar('m');
    nc{ZCOO}.positive = ncchar('down');
    nc{ZCOO}.valid_min =ncshort(0);
    nc{ZCOO}.valid_max =ncshort(2000);
else
    nc{ZCOO} = ncfloat(ZCOO);
    nc{ZCOO}.units =ncchar('kg/m**3');
    nc{ZCOO }.positive = ncchar('down');
    nc{ZCOO}.valid_min =ncshort(0);
    nc{ZCOO}.valid_max =ncshort(60);

end


FILLVALUE        = 32767;

switch PARSAV(1:4)

    case {'TEMP'}
        LONG_NAME = 'Temperature';
        UNITS = 'degree_Celsius';
        VALID_MIN = -3;
        VALID_MAX = 40;
        ADD_OFFSET = 20;
        SCALE_FACT = 0.001d0;

        ER_LONG_NAME = 'Error on temperature (percent variance)';
        ER_UNITS = 'percent of a priori variance';
        ER_VALID_MIN = 0;
        ER_VALID_MAX = 100;

    case {'PSAL'}
        LONG_NAME = 'Salinity';
        UNITS = 'psu';
        VALID_MIN = 0;
        VALID_MAX = 60;
        ADD_OFFSET = 30;
        SCALE_FACT = 0.001d0;

        ER_LONG_NAME = 'Error on salinity (percent variance)';
        ER_UNITS = 'percent of a priori variance';
        ER_VALID_MIN = 0;
        ER_VALID_MAX = 100;

    case {'TPOT'}
        LONG_NAME = 'Potential Temperature';
        UNITS = 'degree_Celsius';
        VALID_MIN = -3;
        VALID_MAX = 40;
        ADD_OFFSET = 20;
        SCALE_FACT = 0.001d0;


    case {'SIG0'}
        LONG_NAME = 'Potential density anomaly';
        UNITS = 'kg/m**3';
        VALID_MIN = 0;
        VALID_MAX = 60;
        ADD_OFFSET = 20;
        SCALE_FACT = 0.001d0;

    case {'SIGI'}
        LONG_NAME = 'In-Situ density anomaly';
        UNITS = 'kg/m**3';
        VALID_MIN = 0;
        VALID_MAX = 60;
        ADD_OFFSET = 20;
        SCALE_FACT = 0.001d0;

    case {'DEPH'}
        LONG_NAME = 'Depth of Potential density anomaly';
        UNITS = 'm';
        VALID_MIN = 0;
        VALID_MAX = 2000;
        ADD_OFFSET = 0;
        SCALE_FACT = 1.00d0;
	
    case {'PRES'}
        LONG_NAME = 'Pressure';
        UNITS = 'decibars';
        VALID_MIN = 0;
        VALID_MAX = 2100;
        ADD_OFFSET = 0;
        SCALE_FACT = 1.00d0;
end

nc{PARSAV} = ncshort('time', ZCOO,'latitude','longitude');   
nc{PARSAV}.long_name    = ncchar(LONG_NAME );  
nc{PARSAV}.units        = ncchar(UNITS);
nc{PARSAV}.valid_min    = ncfloat(VALID_MIN);
nc{PARSAV}.valid_max    = ncfloat(VALID_MAX);
nc{PARSAV}.FillValue_   = ncfloat(FILLVALUE);
nc{PARSAV}.add_offset   = ncfloat(ADD_OFFSET);
nc{PARSAV}.scale_factor = ncfloat(SCALE_FACT);
nc{PARSAV}.comment = ncchar('Estimated by optimal interpolation');


%   Writes coordinates
%  -------------------
nc{'latitude'}(:)  = latitude;
nc{'longitude'}(:) = longitude;
nc{ZCOO}(:)     = depth;
nc{'time'}(:)   = jul_rel;


% Writes parameter
% ----------------
if ~isempty(FIELD)
    inok=find(~finite(FIELD));
    FIELD = FIELD - ADD_OFFSET;
    FIELD = FIELD/SCALE_FACT;
    if ~isempty(inok)
        FIELD(inok)   = FILLVALUE;
    end
    nc{PARSAV}(:,:,:,:) = FIELD;
end

endef(nc);
close(nc);






