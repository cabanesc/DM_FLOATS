% -========================================================
%   USAGE :   verif_prof1(floatname,dacname,numcycle,config_campaign_file,varargin)
%   PURPOSE : compare first descending or ascending profile to ctd reference made at float launch
% -----------------------------------
%   INPUT :
%    floatname  (char)  e.g. '690258'
%    dacname    (char) e.g.  'coriolis'
%    numcycle   (float array)  e.g. 1 (first cycle 001) see option DIRECTION to consider descending profiles
%    config_campaign_file (char)  config file for the campaign e.g 'config_bocats.txt')
%
%   OPTIONNAL INPUT : 
%    'DIRECTION'     (char)    'D' (default): or 'A' to indicate the direction of the profile you want to plot
%    'DATATYPE'      (char)  'raw' (default): if raw PARAM are used; 'adj' if PARAM_ADJ are used  
%    'NB_PROF'       (float)  1    (default): number of profiles before and after "numcycle"  used for the
%                                             plots    
%     'DATAREP'   (char)   default 'DIR_FTP'  directory where nc files are looked for  or 'DIR_DM_FILES'
%    used for interpolating ctd data on float levels:
%    'LIM_SURF_DEEP' (float)  500  (default): Limit depth (m) between surface and deep layers
%    'STEP_SURF'     (float)   10  (default): Distance (m) between measurements in the surface layer
%    'STEP_DEEP'     (float)   20  (default): Distance (m) between measurements in the deep  layer
%  
%    used for plot
%    'ZOOM'           (float) 1000 (default): zoom on  layers below 1000db
% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  :  created: N. David (ENSIETA)- AOUT 2007
%               revised: C. Cabanes - 2016
%   CALLED SUBROUTINES: routine_profil1_m
% ========================================================
function VERIF_PROFIL1(floatname,dacname,numcycle,config_campaign_file,varargin)
close all

% Commented by T. Reynaud: 11/09/2020
%init_path

if iscell(floatname)==0
floatname=cellstr(floatname);
end
if iscell(dacname)==0
dacname=cellstr(dacname);
end
if iscell(numcycle)==0
numcycle={numcycle};
end

if length(floatname)~=length(dacname)|length(floatname)~=length(numcycle)
    error('floatname/dacname and or numcycle do not have the same length')
end


CONFIG=load_configuration('config.txt');

CAMPAIGN=load_configuration(config_campaign_file);

n=length(varargin);

if n/2~=floor(n/2)
    error('check the imput arguments')
end

f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);

% default CONFIG
PARAM.DIRECTION='D';
PARAM.LIM_SURF_DEEP=500;
PARAM.STEP_SURF=10;
PARAM.STEP_DEEP=20;
PARAM.DATATYPE='raw'; 
PARAM.NB_PROF=1;
PARAM.ZOOM=1000;
PARAM.DATAREP='DIR_FTP';

% Input CONFIG
if isfield(s,'DATATYPE')==1;PARAM.DATATYPE=s.DATATYPE;end;
if isfield(s,'LIM_SURF_DEEP')==1;PARAM.LIM_SURF_DEEP=s.LIM_SURF_DEEP;end;
if isfield(s,'STEP_SURF')==1;PARAM.STEP_SURF=s.STEP_SURF;end;
if isfield(s,'STEP_DEEP')==1;PARAM.STEP_DEEP=s.STEP_DEEP;end;
if isfield(s,'DIRECTION')==1;PARAM.DIRECTION=s.DIRECTION;end;
if isfield(s,'NB_PROF')==1;PARAM.NB_PROF=s.NB_PROF;end;
if isfield(s,'ZOOM')==1;PARAM.ZOOM=s.ZOOM;end;
if isfield(s,'DATAREP')==1;PARAM.DATAREP=s.DATAREP;end;
%--------------------------------------------------------------------------
% Entree des parametres
%--------------------------------------------------------------------------

%if exist('auto','var') == 0 % Cas d'execution automatique via routine
%      camp_name = input('Quel est le nom de la campagne oceanographique ? ','s');
%      campagne = input('Quel est le nom du fichier .mat de la campagne (ovid04/ovid06/...) ?','s');  
%      float = input('Quel(s) flotteur(s) ? ');        
%      palier = input('Quel est le palier ? ');
%      pas_surface = input('Quel est le pas en surface ? ');
%      pas_fond = input('Quel est le pas en-dessous du palier ? ');
    

    float = 6901720;
    palier = 500;
    pas_surface = 10;
    pas_fond = 20;
    
    
%end

% 6901631, données <100m
% 6901632, données <0m 
%--------------------------------------------------------------------------
% Traitement
%--------------------------------------------------------------------------

for ifloat = 1:length(floatname)
    flt_name = strtrim(floatname{ifloat});;
	dacname_i =strtrim(dacname{ifloat});
    tab = numcycle{ifloat};
    disp([datestr(now) ' Working on float ' flt_name])
    routine_profil1_m(flt_name,dacname_i,tab,PARAM,CONFIG,CAMPAIGN);
end

%clear all
%close all
