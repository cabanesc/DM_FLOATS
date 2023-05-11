% -========================================================
%   USAGE :   RUN_CPCOR(floatname,dacname,numcycle,config_campaign_file,varargin)
%   PURPOSE : compare first descending or ascending profile to ctd reference made at float launch
% -----------------------------------
%   INPUT :
%    floatname  (char)  e.g. '690258'
%    dacname    (char) e.g.  'coriolis'
%    numcycle   (float array)  e.g. 1 (first cycle 001) see option DIRECTION to consider descending profiles
%    config_campaign_file (char)  config file for the campaign e.g 'config_bocats.txt')
%
%   OPTIONNAL INPUT : 
%    'DIRECTION'     (char)    'A' (default): or 'D' to indicate the direction of the profile you want to plot 
%    used for interpolating ctd data on float levels:
%    'MIN_DEPTH' (float)  1000  (default): Limit depth (m) to be considered for optimization
%    'CORRECT_MINPRES' (float)   2 (default) : minimum pressure difference between two consecutive cycles for which the pressure is corrected by shifting the surface pressure by 1 cycle

% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  : created C. Cabanes - 2020
% ========================================================
function RUN_CPCOR(floatname,dacname,numcycle,config_campaign_file,varargin)
%close all

%init_path

CONFIG=load_configuration('config.txt');

CAMPAIGN=load_configuration(config_campaign_file);
CAMPAIGN.PROF1_DIRECTORY = [CONFIG.DIR_PLOT 'verif_profil1/'];
n=length(varargin);

if n/2~=floor(n/2)
    error('check the imput arguments')
end

f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);

PARAM.MIN_DEPTH=1000;
PARAM.DIRECTION='D';
PARAM.CORRECT_MINPRES=2;

if isfield(s,'MIN_DEPTH')==1;PARAM.MIN_DEPTH=s.MIN_DEPTH;end;
if isfield(s,'DIRECTION')==1;PARAM.DIRECTION=s.DIRECTION;end;
if isfield(s,'CORRECT_MINPRES')==1;PARAM.CORRECT_MINPRES=s.CORRECT_MINPRES;end;

CAMPAIGN.DACNAME=dacname;
CAMPAIGN.FLOAT_SOURCE_NETCDF=CONFIG.DIR_FTP;
CAMPAIGN.DATA_DIRECTORY=CONFIG.DIR_CAMPAIGN;
flt_name = floatname;

fw=fopen([CONFIG.DIR_PLOT 'verif_profil1/'  flt_name '/cpcor_table.csv'],'w')

disp([datestr(now) ' Working on float ' flt_name])
% lo_system_configuration.ZOOM=PARAM.MIN_DEPTH;
%lo_system_configuration.ZOOM='2000';
%--------------------------------------------------------------------------
% Préparation du répertoire de sauvegarde
%--------------------------------------------------------------------------
dir_enregistre = [CONFIG.DIR_PLOT 'verif_profil1/'  flt_name];
if ~exist(dir_enregistre,'dir')
	mkdir(dir_enregistre);
end

routine_profil1_cpcor_prescor(fw,numcycle,flt_name,CAMPAIGN,PARAM);



%clear all
%close all
