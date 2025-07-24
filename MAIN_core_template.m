%    MAIN program to process a CORE float in DM
%    see MAIN_deep_template for deep argo floats.
%
%    Plots several diagnostic figures (theta/S, PSAL, TEMP, SIG0 sections, float trajectory).
%    GUI to view Argo profiles and interactively change the value of the quality flags (corrected flags are reported in netcdf files)
%    Compares an Argo profile with the closest reference profiles.
%    Compares an Argo profile with the reference profile made at the float launch.
%    Analyzes and calculates the salinity correction (using OWC software)
%    Writes DM files with adjusted values. Files are ready to be submitted to the DAC.
%    Verify the adjusted values
%    Produce a report template
%
%    Follow instructions How to install? in readme.md
%    Modify config.txt and rdir below. modify init_path.m if necessary.
%    
close all;
clear all;
curr_dir=pwd;
curr_path=path;

% full path of the DM_FLOATS directory: {REP_CODES}/DM_FLOATS/
rdir='/home/lops/users/ccabanes/dvlpRD/Argo/TD/git/DM_FLOATS/'; % CC
%rdir='/Users/thierry_reynaud/IFREMER/MATLAB/DM_FLOATS_TR/';        % TR
eval(fullfile('cd ',rdir));

% Float WMO numbers for example
%1901210/	3902122/	6901601/	6902803/	6902811/	6902882/
%3900515/	6900807/	6901603/	6902808/	6902818/	6903246/
%3900516/	6901004/	6901758/	6902810/	6902881/	6903249/

floatname = '6902753';
dacname = 'coriolis';
numconfig_ow = 149;                  
% available config_{numconfig_ow}.txt files are in ./LPO_CODES_ATLN_NEW/ow_config/
% 149 : Classical North Atlantic config, ARGO reference database is used
% 129 : Classical North Atlantic config, ARGO  and CTD reference databases are used
%  39 : Classical North Atlantic config, CTD reference databases is used
config_campaign=''; % campaign config files are in VERIF_PROF1/config_campagne/
% available config_campaign files:
%config_bocats.txt    config_geovide.txt  config_ovide18.txt  config_rrex15.txt
%config_catarina.txt  config_ovide10.txt   config_pedro.txt    config_rrex17.txt
%config_oblady18.txt  config_m16420.txt  config_msm9420.txt  config_bocats21.txt
%config_bocats23.txt
%config_arcticgo21.txt config_arcticgo22.txt config_arcticgo23.txt
name_campaign  = 'OVIDE';          % title of the DM report.



n_prof=1; % profile number that we want to compare to the closest profiles in the reference database (see VERIF_FLAG)


irun.LOAD_float      = 1;  % Load data  : copy the file from DIR_FTP_CORIOLIS to DIR_FTP
irun.CORRECT_position= 0;  % Correct the interpolated position (floats in region covered by ice )
irun.PLOTDATA_raw    = 0;  % Preliminary diagnostic plots (theta/S, sections, bathy, flags...)
irun.CORRECT_float   = 0;  % Visualization & correction of flags in netcdf files
irun.VERIF_FLAG      = 0;  % Comparison of a raw Argo profile (n_prof) to the closest profiles in the reference database
irun.VERIF_PROF1_raw = 0;  % Comparison of a raw Argo profile to the launch CTD data
irun.FIND_CLOSE_float= 0;  % Optional: comparison to neighboring Argo profiles (DOI Argo)
irun.OW              = 1;  % OWC correction 
irun.ANALYSE_CAL     = 0;  % Optional: gives information about the OW correction (offset value, slope)
irun.CORRECTIONS     = 0;  % writing D files with OWC correction and corrected flags
irun.PLOTDATA_adj    = 0;  % Diagnostic plots (theta/S, sections, ...) for adjusted data
irun.VERIF_PROF1_adj = 0;  % Comparison of an adjusted Argo profile to the launch CTD data
irun.DOC             = 0;  % Creating the repport

%===========================%
% Load data  : copy the file from DIR_FTP_CORIOLIS to DIR_FTP
if irun.LOAD_float
    rep='CHANGE_FLAGS';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    % LOAD_float(floatname,dacname,'ERASE',1,'ASK',1)  % all options, default values
    LOAD_float(floatname,dacname,'ERASE',1)  % comment this line once the falgs are modified
    eval('cd ..');
    init_path('clear',rep,rdir);
end

%===========================%
% Correction of the interpolated positions (region covered by ice)
if irun.CORRECT_position
    rep='PROG_QC2015';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    %CORR_POSITION(floatname,dacname,'METHOD','original','CORR_NETCDF',0)
    CORR_POSITION(floatname,dacname,'METHOD','sphere','CORR_NETCDF',0)  % to correct in netcdf files: CORR_NETCDF=1 
    eval('cd ..'); 
    init_path('clear',rep,rdir);
end

%===========================%
% Preliminary diagnostic plots (theta/S, sections, bathy,flags...)
if irun.PLOTDATA_raw
    rep='PROG_QC2015';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    % PLOTDATA_func(floatname,dacname,rdir,'DATATYPE','raw','USEFLAG',0,'DATAREP','DIR_FTP','PRINT',1,'PLOT_INVDENS',0,'MAKEPAUSE',0,)    % all options, default values
    PLOTDATA_func(floatname,dacname,rdir,'DATATYPE','raw','PRINT',1,'PLOT_INVDENS',0)   
    PLOT_all_flags(floatname,dacname)
    eval('cd ..');
    init_path('clear',rep,rdir);
end

%===========================%
% Visualization & correction of flags in netcdf files
if irun.CORRECT_float
    rep='CHANGE_FLAGS';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
   %CORRECT_float_flag(floatname,dacname,'VPN',0,'KEEPZOOM',0,'FLAG',1,'NB_PROF',1) % all options, default values
    CORRECT_float_flag(floatname,dacname,'VPN',0,'KEEPZOOM',0,'FLAG',1)% 
    eval('cd ..');
    init_path('clear',rep,rdir);
end

%===========================%
% Comparison of a raw Argo profile (n_prof) to the closest profiles in the reference database
if irun.VERIF_FLAG
    rep='VERIF_FLAG';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    %verif_flag_1profil(floatname,dacname,n_prof,'CONFIG_OW',numconfig_ow','DIRECTION','A','REFERENCE','ctd','REXTEND',50, 'DATATYPE','raw','DATAREP','DIR_FTP','DEPTH_ZOOM',1000,'CHECK_REF',0);  % all options, default values
    verif_flag_1profil(floatname,dacname,n_prof,DEPTH_ZOOM',1000, 'REFERENCE','argo','REXTEND',50);
    eval('cd ..');
    init_path('clear',rep,rdir);
end

%===========================%
% Comparison of a raw Argo profile to the launch CTD data
if irun.VERIF_PROF1_raw
    rep='VERIF_PROF1';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    %VERIF_PROFIL1(floatname,dacname,1,config_campaign,'DIRECTION','D','DATATYPE','raw','NB_PROF',1,'DATAREP','DIR_FTP','NO_FLAG',0,'ZOOM',1000,'CORR_UP_PRES',0);   % all options, default values
    VERIF_PROFIL1(floatname,dacname,1,config_campaign,'NB_PROF',10,'DIRECTION','A','ZOOM',1000);  % Here we compare profile 1A with the launch CTD profile, we plot the 10 first Argo profiles and the upper limit of the vertical axis is 1000 db
    eval('cd ..');
    init_path('clear',rep,rdir);
end

%===========================%
%  OPTIONAL: comparison to neighboring Argo profiles (DOI Argo)
if irun.FIND_CLOSE_float
    rep='COMPARE_FLOAT_REF_TR';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
     % plot a map (at given theta levels) and time series of all observations available in the region (GDAC= DOI ARGO, REF = reference database)
     %PLOT_GDAC_and_FLOAT(floatname,dacname,'UPDATE',0,'PRINT',1,'DATATYPE','raw','TPOT_MIN',4,'TPOT_MAX',4.1,'DEPTH_MIN',0,'VEC_REG',[-30 -20 45 57])
     %PLOT_REF_and_FLOAT(floatname,dacname,'THE_BASES',{'CTD','ARGO'},'UPDATE',0,'PRINT',1,'DATATYPE','raw','TPOT_MIN',4,'TPOT_MAX',4.1,'CYCLE_MIN',1,'CYCLE_MAX',130,'VEC_REG',[-24 -5 22 43],'DEPTH_MIN',0)
     %PLOT_REF_and_FLOAT(floatname,dacname,'THE_BASES',{'CTD','ARGO'},'UPDATE',0,'PRINT',1,'DATATYPE','raw','TPOT_MIN',4,'TPOT_MAX',4.1,'CYCLE_MIN',100,'DEPTH_MIN',0)
     %comparison profile by profile
     %FIND_close_floats(floatname,dacname,'ECARLON',0.2,'ECARLAT',0.1,'ECARDAY',365,'PRES_MIN',1000','UPDATE',1,'PRINT',1,'DATATYPE','raw') % all options, default values
     FIND_close_floats(floatname,dacname,'UPDATE',0,'PRINT',1,'DATATYPE','raw')
    eval('cd ..');
    init_path('clear',rep,rdir);
end

%===========================%
% OW ANALYSE
if irun.OW
    rep='LPO_CODES_ATLN_NEW';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    %MAIN_dmqc_speed(floatname,dacname,numconfig_ow,'RECREATE_MAT','1,'FORCE','','OPTIM',1,'USE_QC',1,'ERASE_MAP',0,'MAKE_PLOT',1,'PLOT_PREVDM',0); % all options, default values
    MAIN_dmqc_speed(floatname,dacname,numconfig_ow,'RECREATE_MAT',1,'USE_QC',1,'ERASE_MAP',0);
    eval('cd ..');
    init_path('clear',rep,rdir);   
    % intermediate repport  and config file for the next steps
    irun.DOC=1;
end

%===========================%
% OPTIONAL: gives information about the OW correction (offset value, slope)
if irun.ANALYSE_CAL
    rep='CORRECTIONS';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    ANALYSE_CAL(floatname,dacname) % A porter TR
    eval('cd ..');
    init_path('clear',rep,rdir);
end

%===========================%
% Writing D files with correction
if irun.CORRECTIONS
    rep='CORRECTIONS';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    MAIN_write_dmqc_files(floatname,dacname) % A porter TR
    eval('cd ..');
    init_path('clear',rep,rdir);
end

%===========================%
% Comparison of an adjusted Argo profile to the launch CTD data
if irun.VERIF_PROF1_adj
    rep='VERIF_PROF1';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    %VERIF_PROFIL1(floatname,dacname,1,config_campaign,'DIRECTION','D','DATATYPE','raw','NB_PROF',1,'DATAREP','DIR_FTP','NO_FLAG',0,'ZOOM',1000,'CORR_UP_PRES',0);   % all options, default values
    VERIF_PROFIL1(floatname,dacname,1,config_campaign,'DIRECTION','A','NB_PROF',10,'DATATYPE','adj','DATAREP','DIR_DM_FILES','ZOOM',1000)
    eval('cd ..');
    init_path('clear',rep,rdir);
end

%===========================%
% Dagnostic plots on adjusted data(theta/S, sections, bathy,flags...)
if irun.PLOTDATA_adj
    rep='PROG_QC2015';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    % PLOTDATA_func(floatname,dacname,rdir,'DATATYPE','raw','USEFLAG',0,'DATAREP','DIR_FTP','PRINT',1,'PLOT_INVDENS',0,'MAKEPAUSE',0,)    % all options, default values
    PLOTDATA_func(floatname,dacname,rdir,'DATATYPE','adj','DATAREP','DIR_DM_FILES','USEFLAG',1,'PRINT',1,'PLOT_INVDENS',0)   
    eval('cd ..');
    init_path('clear',rep,rdir);
end

%===========================%
% FINAL REPORT
if irun.DOC
    rep='DOC';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    %n_prof=0;% A revoir avec CC
    %generate_doc_overleaf(floatname,dacname,numconfig_ow,'SUBTITLE',['Float ' floatname ' (' name_campaign ')'],'PROFREF',n_prof,'PROFFLAG','1D')%    
    generate_doc_overleaf(floatname,dacname,numconfig_ow,'SUBTITLE',['Float ' floatname ' (' name_campaign ')'],'PROFREF',n_prof,'COMP_GDAC',1)
    %generate_doc_overleaf(floatname,dacname,numconfig_ow,'SUBTITLE',['Float ' floatname ' (' name_campaign ')'],'PROFREF',{[40,56,70,86]})
    %generate_doc_overleaf(floatname,dacname,numconfig_ow,'SUBTITLE',['Float ' floatname ' (' name_campaign ')'],'PROFREF',n_prof)
    eval('cd ..');
    init_path('clear',rep,rdir);
    % upload files in ./DOC/OVERLEAF/ on
    % https://www.overleaf.com/ et compile  .tex file
end
%===========================%
% save the current run
if exist ('./TEMPLATES')==0
    mkdir ./TEMPLATES
end

if ~(contains(curr_dir,'TEMPLATES'))&~(contains(curr_path,'TEMPLATES'))
copyfile ('MAIN_core.m', ['./TEMPLATES/MAIN_' floatname '.m'])
else
eval(['cd ' curr_dir])
end

SAUV_RAPP=0; %sauvegarde du rapport
if SAUV_RAPP==1
eval(['!cp ./DOC/OVERLEAF/ebauche_rapport.tex ./DOC/OVERLEAF/rapport_' floatname '.tex'])
eval(['!cp -R  ./DOC/OVERLEAF/rapport_' floatname '.pdf /home5/pharos/argo/DMARGO/data/DM_FILES/Tech_Reports_PSAL/' floatname '.pdf'])
eval(['!chmod 755 /home5/pharos/argo/DMARGO/data/DM_FILES/Tech_Reports_PSAL/' floatname '.pdf'])
% sauvegarde des données
eval(['!rm -rf /home5/pharos/argo/DMARGO/data/DM_FILES/coriolis/' floatname])
eval(['!cp -rf /export/home1/ccabanes/data/DM_FILES/coriolis/' floatname ' /home5/pharos/argo/DMARGO/data/DM_FILES/coriolis/.'])
eval(['!chmod -R 755 /home5/pharos/argo/DMARGO/data/DM_FILES/coriolis/' floatname])
%eval(['!cp -f /export/home1/ccabanes/data/FTP_ARGO/coriolis/' floatname '/profiles/B*.nc ' '/home5/pharos/argo/DMARGO/data/DM_FILES/coriolis/' floatname '/profiles/.'])
eval(['!chmod -R 755 /home5/pharos/argo/DMARGO/data/DM_FILES/coriolis/' floatname])
end
