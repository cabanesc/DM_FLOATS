
% -========================================================
%   USAGE : MAIN_dmqc_speed(tabfloat,tabdac,configow,varargin)
%   PURPOSE : run OWC software for a float or a group of floats, using a specified OW config
% -----------------------------------
%   INPUT :
%    tabfloat  (char or cell of chars -size n_floatsx1)    e.g. '6900258' or {'6900258', '3901954'}
%    tabdac    (char or cell of chars -size n_floatsx1)    e.g. 'coriolis' or {'coriolis', 'bodc'}
%    configow  (floats or cell of floats -size n_floatsx1) e.g.  149       or {149,149}      % config number ow
%
%   OPTIONNAL INPUT :
%    'RECREATE_MAT' (logical)  RECREATE_MAT=1 creates again the source mat file. If RECREATE_MAT=0 (default), uses the existing source mat file.
%     'FORCE'        (char)   'raw' : force create_float_source.m to load PRES TEMP and PSAL
%                             'adjusted': force create_float_source.m to load PRES_ADJUSTED, PSAL_ADJUSTED, TEMP_ADJUSTED
%                             '' (default is empty) => variables loaded are raw PRES, PSAL and TEMP.
%                                                      If PRES is adjusted, variables loaded are PRES_ADJUSTED, raw PSAL calibrated in pressure and raw TEMP.
%    'OPTIM'        (logical)  OPTIM=1 uses optimization toolbox (default), otherwise use LMA function 
%    'USE_PRES_LT'      (float)      e.g 1500 (default is empty)
%    'USE_PRES_GT'      (float)      e.g 2000 (default is empty)
%    'POSTFIX'            (char)      'eg _1500_2000' used to document the set_calseries options (default is empty)
%    'USE_QC'         (logical)   USE_QC=1 (default) uses PSAL_QC and remove PSAL data with QC=4
%                                 USE_QC=0           do not use PSAL_QC (do not remove data with QC=4)  
%     'ERASE_MAP'      (logical)   ERASE_MAP=0 (default) keep the previous  map file ; ERASE_MAP=1 erase the previous map file
%     'MAKE_PLOT'        (logical)   MAKE_PLOT=1 (default) make all plots; MAKE_PLOT=0 does not make plots
% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  : created (2016) ccabanes
%
%   CALLED SUBROUTINES: 
% -------------------------------------
% 
% ========================================================
function MAIN_dmqc_speed(tabfloat,tabdac,configow,varargin)

%init_path

C = load_configuration('config.txt');

% dir input /output
DIR_DATA=C.DIR_DATA;
% Added by T. Reynaud 21.09.2020:
DIR_OWC=C.DIR_OWC;
if exist([DIR_DATA  'float_calib'])==0; mkdir([DIR_DATA  'float_calib']); end;
if exist([DIR_DATA  'float_mapped'])==0; mkdir([DIR_DATA  'float_mapped']); end;
if exist([DIR_DATA  'float_plots'])==0; mkdir([DIR_DATA  'float_plots']); end;
if exist([DIR_DATA  'float_source'])==0; mkdir([DIR_DATA  'float_source']); end;


if iscell(tabfloat)==0;tabfloat=cellstr(tabfloat);end
if iscell(tabdac)==0;tabdac=cellstr(tabdac);end
if iscell(configow)==0;configow={configow};end
if length(tabfloat)>1&length(tabdac)==1
tabdac=repmat(tabdac,1,length(tabfloat)); 
end
if length(tabfloat)>1&length(configow)==1
configow=repmat(configow,1,length(configow)); 
end

n=length(varargin);

% INPUT PARAMETERS
if n/2~=floor(n/2)
    error('check the imput arguments')
end

f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);

% default 
RECREATE_MAT=0;
OPTIM=1;
POSTFIX='';
FORCE='';
USE_PRES_LT=[];
USE_PRES_GT=[];
USE_THETA_LT=[];
USE_THETA_GT=[];
USE_QC=1;
ERASE_MAP=0;
MAKE_PLOT=1;
if isfield(s,'RECREATE_MAT')==1;RECREATE_MAT=s.RECREATE_MAT;end;
if isfield(s,'OPTIM')==1;OPTIM=s.OPTIM;end;
if isfield(s,'POSTFIX')==1;POSTFIX=s.POSTFIX;end;
if isfield(s,'FORCE')==1;FORCE=s.FORCE;end;
if isfield(s,'USE_PRES_GT')==1;USE_PRES_GT=s.USE_PRES_GT;end;
if isfield(s,'USE_PRES_LT')==1;USE_PRES_LT=s.USE_PRES_LT;end;
if isfield(s,'USE_QC')==1;USE_QC=s.USE_QC;end;
if isfield(s,'ERASE_MAP')==1;ERASE_MAP=s.ERASE_MAP;end;
if isfield(s,'MAKE_PLOT')==1;MAKE_PLOT=s.MAKE_PLOT;end;


RECREATE_MAT
% VERSION of the code/reference data:
VERSION_OW = C.VERSION_OWC;


if isempty(POSTFIX)==0
	[a,r]=strtok(fliplr(POSTFIX),'_'); 
	USE_PRES_LT=str2num(fliplr(a));
	[a,r]=strtok((r),'_'); 
	USE_PRES_GT=str2num(fliplr(a));
	if isempty(USE_PRES_GT)
	USE_PRES_GT=[];
	end
	if isempty(USE_PRES_LT)
	USE_PRES_LT=[];
    end
    
    [a,r]=strtok(fliplr(POSTFIX),'T'); 
	USE_THETA_LT=str2num(fliplr(a));
	[a,r]=strtok((r),'T'); 
	USE_THETA_GT=str2num(fliplr(a));
	if isempty(USE_THETA_GT)
	USE_THETA_GT=[];
	end
	if isempty(USE_THETA_LT)
	USE_THETA_LT=[];
	end
end

%PATHS
if OPTIM==0
	if exist('lsqnonlin')
		[r,optim_path]=strtok(fliplr(which('lsqnonlin')),'/');
		 optim_path = fliplr(optim_path);
		 rmpath(optim_path)
	end
end

%  path codes OW dans cet ordre!
addpath([C.DIR_OWC VERSION_OW '/matlab_codes/']) % codes OW originaux
addpath('util/')
addpath('ow_codes_modif/')  % modif codes OW
addpath('ow_codes_modif/Speed/')  % modif codes OW
% sauve le path
path_save=path;

%fid10=fopen('/tmp/convert.sh','w')

for ifloat=1:length(tabfloat)
    
    flt_name=tabfloat{ifloat};
    flt_dir='';
    flt_dir = deblank(flt_dir);
    flt_name = deblank(flt_name);
    flt_name
	
    for numconfig = configow{ifloat};
        
        %thecal='_6901602';
        % reinitialise le path
        path(path_save)
        
        addpath('ow_codes_modif/For_config129/')
        % if numconfig==14
            % addpath('ow_codes_modif/For_config14/')  % REF ARGO + calcul d'un offset seul
        % end
        
        % if numconfig==3 % REF CTD + calcul d'un offset seul
            % addpath('ow_codes_modif/For_config3/')
        % end
        
        % if numconfig==127  % REF CTD + ARGO, selection des donnees historiques  +/- 2yr, calcul d'un offset seul
            % addpath('ow_codes_modif/For_config127/')
        % end
        
        % if numconfig==129 |numconfig==1291|numconfig==1292|numconfig==1293|numconfig==1295|numconfig==149|numconfig==1491 |numconfig==1492|numconfig==1493|numconfig==1495|numconfig==1494|numconfig==14911% REF CTD + ARGO, selection des donnees historiques  +/- map_age_large (voir config), recalcul des erreurs, calcul d'un offset seul
            % addpath('ow_codes_modif/For_config129/')
        % end
         % if numconfig==128 % test pour birgit
            % addpath('ow_codes_modif/For_config128/')
        % end
        % if numconfig==39|numconfig==392|numconfig==395|numconfig==391|numconfig==393  % REF CTD , selection des donnees historiques  +/- map_age_large (voir config), recalcul des erreurs, calcul d'un offset seul
            % addpath('ow_codes_modif/For_config39/')
        % end
        % if numconfig==371  % REF CTD, selection des donnees historiques  +/- 2yr, echelle spatiales (LONG et LAT) divisees par 4
            % addpath('ow_codes_modif/For_config37/')
        % end
        
        thecal='';
        thecal=['_' flt_name];
        if exist(['set_calseries' thecal '.m'])==0
            thecal='';
        end
        
        % charge la config
        % DIW_OWC added by T. Reynaud 21.09.2020
        lo_system_configuration = load_configuration( ['./ow_config/ow_config_' num2str(numconfig) '.txt'] ,flt_name, DIR_DATA, VERSION_OW, DIR_OWC);
        lo_system_configuration.use_pres_gt=USE_PRES_GT;
		lo_system_configuration.use_pres_lt=USE_PRES_LT;
		lo_system_configuration.use_theta_gt=USE_THETA_GT;
		lo_system_configuration.use_theta_lt=USE_THETA_LT;
		
        % creation du repertoire des plots
        dirplot_config = [ DIR_DATA '/float_plots/CONFIG' num2str(numconfig) '/'];
        dirplot = [dirplot_config  flt_name '/'];
        
        if exist(dirplot_config) == 0;   mkdir(dirplot_config);     end;
        if exist(dirplot) == 0;   mkdir(dirplot);            end;
        
        % creation du repertoire du mapping
        
        dirmap = [ DIR_DATA '/float_mapped/CONFIG' num2str(numconfig) '/'];
        if exist(dirmap) == 0; mkdir(dirmap);    end;
        
        % Ajout Cecile : creation auto du repertoire des calib
        
        dircalib = [DIR_DATA 'float_calib/CONFIG' num2str(numconfig) '/'];
        
        if exist(dircalib) == 0;  mkdir(dircalib);    end;
        
        % Ajout Cecile : test si le fichier .mat existe
        auto=1;
        float = flt_name;
        dacname = deblank(tabdac{ifloat});
        
        dirsource = [DIR_DATA '/float_source/'];
		
		if exist([dirsource strtrim(flt_name) '.mat'],'file')==2
			if RECREATE_MAT==1
			
			   eval(['!/bin/rm -f ' dirsource strtrim(flt_name) '.mat'])
			end
		end
		if isempty(FORCE)
        create_float_source(flt_name,dacname,'useqc',USE_QC)
		else
		create_float_source(flt_name,dacname,'force',FORCE,'useqc',USE_QC);
		end
        
        
        
        % copie des fichier map_ si ils ont deja ete crees pour une autre config (sans changement des parametres)
        
        if numconfig==14 % REF ARGO + calcul d'un offset seul
            % si la config 11 a deja ete faite, on copie le fichier map_
            if exist([DIR_DATA '/float_mapped/CONFIG11/map_' flt_name '.mat'])
                copyfile([DIR_DATA '/float_mapped/CONFIG11/map_' flt_name '.mat'],[DIR_DATA '/float_mapped/CONFIG' num2str(numconfig) '/.']);
            end
        end
        
        if numconfig==3 % REF CTD + calcul d'un offset seul
            % si la config 1 a deja ete faite, on copie le fichier map_
            if exist( [DIR_DATA '/float_mapped/CONFIG1/map_' flt_name '.mat'])
                copyfile([DIR_DATA '/float_mapped/CONFIG1/map_' flt_name '.mat'],[DIR_DATA '/float_mapped/CONFIG' num2str(numconfig) '/.']);
            end
        end
        % plot_corr_deep_fct_pres( flt_dir, flt_name, lo_system_configuration );
        
        disp([datestr(now) ' Working on : flotteur ' flt_name])
        eval(['!/bin/rm -f ' DIR_DATA 'float_calib/CONFIG' num2str(numconfig) '/calseries_' num2str(flt_name) '.mat'])
        texte=['set_calseries' thecal '( flt_dir, flt_name, lo_system_configuration );']
        eval(texte)
        %set_calseries( flt_dir, flt_name, lo_system_configuration ); 
        
        % update_salinity_mapping( flt_dir, flt_name, lo_system_configuration );
        % plus rapide
        if ERASE_MAP==1
        eval(['!/bin/rm -f ' DIR_DATA 'float_mapped/CONFIG' num2str(numconfig) '/map_' num2str(flt_name) '.mat'])
        end
        update_salinity_mapping_speed( flt_dir, flt_name, lo_system_configuration );
        
        calculate_piecewisefit( flt_dir, flt_name, lo_system_configuration );
		
		if isempty(POSTFIX)==0
		   eval(['!/bin/cp -f ' dircalib 'cal_' flt_name '.mat ' dircalib 'cal_' flt_name POSTFIX '.mat'])
		   eval(['!/bin/cp -f ' dircalib 'calseries_' flt_name '.mat ' dircalib 'calseries_' flt_name POSTFIX '.mat'])
        end
        if MAKE_PLOT==1
            plot_diagnostics_ow( flt_dir, flt_name, lo_system_configuration );
            %keyboard
            plot_diagnostics_ow_figure9_4( flt_dir, flt_name, lo_system_configuration,dacname,C);
            
            plotfile=dir([DIR_DATA '/float_plots/CONFIG' num2str(numconfig) '/' flt_name '/*.eps']);
        end
%          for k=1:length(plotfile)
%              thefile=[DIR_DATA '/float_plots/CONFIG' num2str(numconfig) '/' flt_name '/' plotfile(k).name];
%              thestring=['!ps2pdf ' thefile ' ' strrep(thefile,'.eps',[POSTFIX '.pdf'])];
%              eval(thestring)
%              thestring=['!rm -f ' thefile ];
%              eval(thestring)
%          end
         
        
    end % fin boucle sur configs
   
end

