%-------------------------------------------------------------------------------
%c
%c mltnc_argofile_create    - Creation du fichier multistations au format NetCDF 		
%c
%-------------------------------------------------------------------------------
%  version:
%  --------
%  1.01   Creation              			30/04/03  E. Autret
%       
% Modification Carole GRIT / Fevrier 2004 pour Etude SPMW
%
%-------------------------------------------------------------------------------
function [msg_error] = mltnc_argofile_create_pres(ncfile_name,software,n_prof,n_param,nb_nivstd);
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------


% ------------------------
%    Create the new file
% ------------------------
nc = netcdf(ncfile_name, 'clobber');
msg_error = 'ok';
if isempty(nc)
   msg_error = ['mltnc_argofile_create : Problème d''ouverture du fichier : ' ncfile_name];
end



% ----------------------------
%    Defines the dimensions
% ----------------------------
%  Fixed dimensions
nc('DATE_TIME') = 14;  %date_time convention is YYYYMMDDHHMISS
nc('STRING256') = 256;
nc('STRING64')  = 64;
nc('STRING32')  = 32;		
nc('STRING16')  = 16;
nc('STRING8')   = 8;
nc('STRING4')   = 4;
nc('STRING2')   = 2;

%  Variable dimensions
nc('N_PROF')     = n_prof;	 %  Number of profiles
nc('N_PARAM')    = n_param ;  %  Maximum number of parameters
nc('N_LEVELS')   = nb_nivstd;    %  Maximum number of Z levels

% ----------------------------
%    Defines global attributes
% ----------------------------

nc.Last_update   = ncchar(datestr(now, 0));
nc.SoftwareVersion       = ncchar(software);

% ----------------------------------------
%   Variables and attributes:
% ----------------------------------------
%
nc{'DATA_TYPE'} = ncchar('STRING16');%% 16 elements.
nc{'DATA_TYPE'}.comment  = ncchar('Data type') ;

nc{'FORMAT_VERSION'} = ncchar('STRING4');%% 4 elements.
nc{'FORMAT_VERSION'}.comment   =ncchar('File format version') ;

nc{'REFERENCE_DATE_TIME'} = ncchar('DATE_TIME'); %% 14 elements.
nc{'REFERENCE_DATE_TIME'}.comment = ncchar('Date of reference for Julian days');
nc{'REFERENCE_DATE_TIME'}.conventions = ncchar('YYYYMMDDHHMISS');

nc{'PROJECT_NAME'} = ncchar('N_PROF', 'STRING64'); %% 586496 elements.
nc{'PROJECT_NAME'}.comment = ncchar('Name of the project');

nc{'PI_NAME'} = ncchar('N_PROF', 'STRING64'); %% 586496 elements.
nc{'PI_NAME'}.comment = ncchar('Name of the principal investigator');


nc{'PLATFORM_NUMBER'} = ncchar('N_PROF', 'STRING8'); %% 73312 elements.
nc{'PLATFORM_NUMBER'}.long_name = ncchar('Float unique identifier');
nc{'PLATFORM_NUMBER'}.conventions = ncchar('WMO float identifier: QA9IIIII');

nc{'CYCLE_NUMBER'} = nclong('N_PROF'); %% 9164 elements.
nc{'CYCLE_NUMBER'}.long_name = ncchar('Float cycle number');
nc{'CYCLE_NUMBER'}.conventions = ncchar('0..N, 0 : launch cycle (if exists), 1 : first complete cycle');
nc{'CYCLE_NUMBER'}.FillValue_ = nclong(99999);

nc{'DIRECTION'} = ncchar('N_PROF'); %% 9164 elements.
nc{'DIRECTION'}.long_name = ncchar('Direction of the station profiles');
nc{'DIRECTION'}.conventions = ncchar('A: ascending profiles, D: descending profiles');

nc{'DATA_CENTRE'} = ncchar('N_PROF', 'STRING2'); %% 18328 elements.
nc{'DATA_CENTRE'}.long_name = ncchar('Data centre in charge of float data processing');
nc{'DATA_CENTRE'}.conventions = ncchar('GTSPP table');

%nc{'DATE_CREATION'} = ncchar('DATE_TIME'); %% 14 elements.
%nc{'DATE_CREATION'}.comment = ncchar('Date of file creation');
%nc{'DATE_CREATION'}.conventions = ncchar('YYYYMMDDHHMISS');

%nc{'DATE_UPDATE'} = ncchar('DATE_TIME'); %% 14 elements.
%nc{'DATE_UPDATE'}.long_name = ncchar('Date of update of this file');
%nc{'DATE_UPDATE'}.conventions = ncchar('YYYYMMDDHHMISS');

nc{'DC_REFERENCE'} = ncchar('N_PROF', 'STRING32'); %% 146624 elements.
nc{'DC_REFERENCE'}.long_name = ncchar('Station unique identifier in data centre');
nc{'DC_REFERENCE'}.conventions = ncchar('Data centre convention');

nc{'DATA_STATE_INDICATOR'} = ncchar('N_PROF', 'STRING4'); %% 36656 elements.
nc{'DATA_STATE_INDICATOR'}.long_name = ncchar('Degree of processing the data have passed through');
nc{'DATA_STATE_INDICATOR'}.conventions = ncchar('OOPC table');

nc{'DATA_MODE'} = ncchar('N_PROF'); 
nc{'DATA_MODE'}.long_name = ncchar('Delayed mode or real time data');
nc{'DATA_MODE'}.conventions = ncchar('R : real time; D : delayed mode');

nc{'INST_REFERENCE'} = ncchar('N_PROF', 'STRING64'); 
nc{'INST_REFERENCE'}.long_name = ncchar('Instrument type');
nc{'INST_REFERENCE'}.conventions = ncchar('Brand, type, serial number');

nc{'WMO_INST_TYPE'} = ncchar('N_PROF', 'STRING4');
nc{'WMO_INST_TYPE'}.long_name = ncchar('Coded instrument type');
nc{'WMO_INST_TYPE'}.conventions = ncchar('WMO code table 1770 - instrument type');

nc{'JULD'} = ncdouble('N_PROF'); 
nc{'JULD'}.long_name = ncchar('Julian day (UTC) of the station relative to REFERENCE_DATE_TIME');
nc{'JULD'}.units = ncchar('days since 1950-01-01 00:00:00 UTC');
nc{'JULD'}.conventions = ncchar('Relative julian days with decimal part (as parts of day)');
nc{'JULD'}.FillValue_ = ncdouble(999999);

%nc{'JULD_QC'} = ncchar('N_PROF'); 
%nc{'JULD_QC'}.long_name = ncchar('Quality on Date and Time');
%nc{'JULD_QC'}.conventions = ncchar('Q where Q =[0-9]');
%nc{'JULD_QC'}.FillValue_ = ncchar('0');

nc{'LATITUDE'} = ncdouble('N_PROF');
nc{'LATITUDE'}.long_name = ncchar('Latitude of the station, best estimate');
nc{'LATITUDE'}.units = ncchar('degree_north');
nc{'LATITUDE'}.FillValue_ = ncdouble(99999);
nc{'LATITUDE'}.valid_min = ncdouble(-90);
nc{'LATITUDE'}.valid_max = ncdouble(90);

nc{'LONGITUDE'} = ncdouble('N_PROF');
nc{'LONGITUDE'}.long_name = ncchar('Longitude of the station, best estimate');
nc{'LONGITUDE'}.units = ncchar('degree_east');
nc{'LONGITUDE'}.FillValue_ = ncdouble(99999);
nc{'LONGITUDE'}.valid_min = ncdouble(-180);
nc{'LONGITUDE'}.valid_max = ncdouble(180);

%nc{'PRES'} = ncfloat('N_LEVELS'); 

nc{'PRES'} = ncfloat('N_PROF', 'N_LEVELS'); 
nc{'PRES'}.long_name = ncchar('Pressure');
nc{'PRES'}.FillValue_ = ncfloat(99999);
nc{'PRES'}.units = ncchar('decibar');
nc{'PRES'}.valid_min = ncfloat(0);
nc{'PRES'}.valid_max = ncfloat(15000);
nc{'PRES'}.comment = ncchar('In situ measurement');
%nc{'PRES'}.C_format = ncchar('%7.1f');
%nc{'PRES'}.FORTRAN_format = ncchar('F7.1');
%nc{'PRES'}.resolution = ncfloat(0.100000001490116);

nc{'PRES_ERR'} = ncfloat('N_PROF', 'N_LEVELS'); 
nc{'PRES_ERR'}.long_name = ncchar('Error on pressure');
nc{'PRES_ERR'}.FillValue_ = ncfloat(99999);
nc{'PRES_ERR'}.units = ncchar('decibar');
nc{'PRES_ERR'}.valid_min = ncfloat(0);
nc{'PRES_ERR'}.valid_max = ncfloat(0.1);


nc{'PRES_QC'} = ncchar('N_PROF', 'N_LEVELS'); 
nc{'PRES_QC'}.long_name = ncchar('Quality on pressure');
nc{'PRES_QC'}.conventions = ncchar('Q where Q =[0-5]');
nc{'PRES_QC'}.FillValue_ = ncchar('0');


nc{'TEMP'} = ncfloat('N_PROF', 'N_LEVELS'); 
nc{'TEMP'}.long_name = ncchar('Temperature in situ T90 scale (interpolated on PRES levels)');
nc{'TEMP'}.FillValue_ = ncfloat(99999);
nc{'TEMP'}.units = ncchar('degree_Celsius');
nc{'TEMP'}.valid_min = ncfloat(-3);
nc{'TEMP'}.valid_max = ncfloat(40);
nc{'TEMP'}.comment = ncchar('In situ measurement');

nc{'TEMP_ERR'} = ncfloat('N_PROF', 'N_LEVELS'); 
nc{'TEMP_ERR'}.long_name = ncchar('Error on interpolated temperature');
nc{'TEMP_ERR'}.FillValue_ = ncfloat(99999);
nc{'TEMP_ERR'}.units = ncchar('degree_Celsius');
nc{'TEMP_ERR'}.valid_min = ncfloat(0);
nc{'TEMP_ERR'}.valid_max = ncfloat(0.1);


nc{'TEMP_QC'} = ncchar('N_PROF', 'N_LEVELS'); 
nc{'TEMP_QC'}.long_name = ncchar('Quality on interpolated temperature');
nc{'TEMP_QC'}.conventions = ncchar('Q where Q =[0-5]');
nc{'TEMP_QC'}.FillValue_ = ncchar('0');



nc{'PSAL'} = ncfloat('N_PROF', 'N_LEVELS');
nc{'PSAL'}.long_name = ncchar('Practical salinity (interpolated on PRES levels)');
nc{'PSAL'}.FillValue_ = ncfloat(99999);
nc{'PSAL'}.units = ncchar('psu');
nc{'PSAL'}.valid_min = ncfloat(0);
nc{'PSAL'}.valid_max = ncfloat(60);
nc{'PSAL'}.comment = ncchar('In situ measurement');
  

nc{'PSAL_ERR'} = ncfloat('N_PROF', 'N_LEVELS'); 
nc{'PSAL_ERR'}.long_name = ncchar('Error on interpolated practical salinity');
nc{'PSAL_ERR'}.FillValue_ = ncfloat(99999);
nc{'PSAL_ERR'}.units = ncchar('psu');
nc{'PSAL_ERR'}.valid_min = ncfloat(0);
nc{'PSAL_ERR'}.valid_max = ncfloat(0.1);


nc{'PSAL_QC'} = ncchar('N_PROF', 'N_LEVELS'); 
nc{'PSAL_QC'}.long_name = ncchar('Quality on interpolated practical salinity');
nc{'PSAL_QC'}.conventions = ncchar('Q where Q =[0-5]');
nc{'PSAL_QC'}.FillValue_ = ncchar('0');


  
nc{'TPOT'} = ncfloat('N_PROF', 'N_LEVELS'); 
nc{'TPOT'}.long_name = ncchar('Potential temperature  (interpolated on PRES levels)');
nc{'TPOT'}.FillValue_ = ncfloat(99999);
nc{'TPOT'}.units = ncchar('degree_Celsius');
nc{'TPOT'}.valid_min = ncfloat(-3);
nc{'TPOT'}.valid_max = ncfloat(40);
nc{'TPOT'}.comment = ncchar('Deduced from in situ measurement');

nc{'TPOT_ERR'} = ncfloat('N_PROF', 'N_LEVELS'); 
nc{'TPOT_ERR'}.long_name = ncchar('Error on interpolated potential temperature');
nc{'TPOT_ERR'}.FillValue_ = ncfloat(99999);
nc{'TPOT_ERR'}.units = ncchar('degree_Celsius');
nc{'TPOT_ERR'}.valid_min = ncfloat(0);
nc{'TPOT_ERR'}.valid_max = ncfloat(0.1);


nc{'TPOT_QC'} = ncchar('N_PROF', 'N_LEVELS'); 
nc{'TPOT_QC'}.long_name = ncchar('Quality on interpolated potential temperature');
nc{'TPOT_QC'}.conventions = ncchar('Q where Q =[0-5]');
nc{'TPOT_QC'}.FillValue_ = ncchar('0');



nc{'SIG0'} = ncfloat('N_PROF', 'N_LEVELS'); 
nc{'SIG0'}.long_name = ncchar('Density Anomaly referenced to P=0 (interpolated on PRES levels)');
nc{'SIG0'}.FillValue_ = ncfloat(99999);
nc{'SIG0'}.units = ncchar('kg/m3');
nc{'SIG0'}.valid_min = ncfloat(0);
nc{'SIG0'}.valid_max = ncfloat(60);
nc{'SIG0'}.comment = ncchar('Deduced from in situ measurement');

nc{'SIG0_ERR'} = ncfloat('N_PROF', 'N_LEVELS'); 
nc{'SIG0_ERR'}.long_name = ncchar('Error on interpolated density anomaly');
nc{'SIG0_ERR'}.FillValue_ = ncfloat(99999);
nc{'SIG0_ERR'}.units = ncchar('kg/m3');
nc{'SIG0_ERR'}.valid_min = ncfloat(0);
nc{'SIG0_ERR'}.valid_max = ncfloat(0.1);


nc{'SIG0_QC'} = ncchar('N_PROF', 'N_LEVELS'); 
nc{'SIG0_QC'}.long_name = ncchar('Quality on interpolated density anomaly');
nc{'SIG0_QC'}.conventions = ncchar('Q where Q =[0-5]');
nc{'SIG0_QC'}.FillValue_ = ncchar('0');



nc{'VORP'} = ncfloat('N_PROF', 'N_LEVELS'); 
nc{'VORP'}.long_name = ncchar('Planetary vorticity (interpolated on PRES levels)');
nc{'VORP'}.FillValue_ = ncfloat(99999);
nc{'VORP'}.units = ncchar('m-1.s-1');
nc{'VORP'}.valid_min = ncfloat(-1.e-12);
nc{'VORP'}.valid_max = ncfloat(+1.e-12);
nc{'VORP'}.comment = ncchar('Deduced from in situ measurement');


nc{'VORP_ERR'} = ncfloat('N_PROF', 'N_LEVELS'); 
nc{'VORP_ERR'}.long_name = ncchar('Error on interpolated potential vorticity');
nc{'VORP_ERR'}.FillValue_ = ncfloat(99999);
nc{'VORP_ERR'}.units = ncchar('m-1.s-1');
nc{'VORP_ERR'}.valid_min = ncfloat(-1.e-12);
nc{'VORP_ERR'}.valid_max = ncfloat(+1.e-12);


nc{'VORP_QC'} = ncchar('N_PROF', 'N_LEVELS'); 
nc{'VORP_QC'}.long_name = ncchar('Quality on interpolated potential vorticity');
nc{'VORP_QC'}.conventions = ncchar('Q where Q =[0-5]');
nc{'VORP_QC'}.FillValue_ = ncchar('0');


endef(nc);
close(nc);
