%  MAIN program to process a float in DM
% chargement des données, plots preliminaires (theta/S, sections, bathy...),visualisation & correction des flags
% comparaison aux profils voisins, methode OWC, creation des fichiers D avec correction, verification des données ajustees. 
%  modify config.txt and rdir below. modify init_path.m if necessary.
close all;
clear all;
% full path of the DM_FLOATS directory: rdir
rdir='/home1/homedir5/perso/ccabanes/dvlpRD/Argo/TD/git/DM_FLOATS/'; % CC
%rdir='/Users/thierry_reynaud/IFREMER/MATLAB/DM_FLOATS_TR/';        % TR
eval(fullfile('cd ',rdir));

% Float WMO numbers for example
%1901210/	3902122/	6901601/	6902803/	6902811/	6902882/
%3900515/	6900807/	6901603/	6902808/	6902818/	6903246/
%3900516/	6901004/	6901758/	6902810/	6902881/	6903249/

floatname = '3901919';
dacname = 'coriolis';
numconfig_ow = 1494;                  % available ow configuration files are in LPO_CODES_ATLN_NEW/ow_config/ :
% 149 : Classical North Atlantic config, ARGO reference database is used
% 129 : Classical North Atlantic config, ARGO  and CTD reference databases are used
%  39 : Classical North Atlantic config, CTD reference databases is used
config_campaign=''; % campaign config files are in VERIF_PROF1/config_campagne/
% available config_campaign files:
%config_bocats.txt    config_geovide.txt  config_ovide18.txt  config_rrex15.txt
%config_catarina.txt  config_ovide10.txt   config_pedro.txt    config_rrex17.txt
name_campaign  = 'MOCCA ';          % titre pour le rapport



n_prof=10; % numero de profil que l'on veut verifier avec VERIF_FLAG


irun.LOAD_float      = 0;  % chargement des données: copie les fichiers netcdf depuis DIR_FTP_CORIOLIS vers DIR_FTP
irun.PLOTDATA_raw    = 0;  % plots preliminaires (theta/S, sections, bathy...)
irun.CORRECT_float   = 0;  % visualisation & correction des flags
irun.VERIF_FLAG      = 0;  % comparaison d'un profil Argo (n_prof) aux profils les plus proches de la base de reference
irun.VERIF_PROF1_raw = 0;  % comparaison au profil CTD de mise a l'eau
irun.FIND_CLOSE_float= 1;  % optionnel : comparaison aux profils Argo voisins 
irun.OW              = 0;  % calcul correction OWC
irun.CORRECTIONS     = 0;  % ecriture des fichiers  D  avec correction
irun.PLOTDATA_adj    = 0;  % verification des donnees ajustees (theta/S, sections)
irun.VERIF_PROF1_adj = 0;  % verification des donnees ajustees (comparaison au profil CTD de mise a l'eau)
irun.DOC             = 0;  % creation du rapport

%===========================%
% chargement des données depuis le ftp
if irun.LOAD_float
    rep='CHANGE_FLAGS';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    LOAD_float(floatname,dacname,'ERASE',1)  % commenter par securite une fois que les flags sont modifies
    eval('cd ..');
    init_path('clear',rep,rdir);
    % premier diag, en particulier inversion de densité
end
%keyboard
if irun.PLOTDATA_raw
    rep='PROG_QC2015';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    PLOTDATA_func(floatname,dacname,rdir,'DATATYPE','raw','PRINT',1)   % attention matlab en display
    eval('cd ..');
    init_path('clear',rep,rdir);
end


%===========================%
% modifocation des flags si necessaire
if irun.CORRECT_float
    rep='CHANGE_FLAGS';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    CORRECT_float_flag(floatname,dacname,'VPN',0)% RAS
    
    eval('cd ..');
    init_path('clear',rep,rdir);
end
%===========================%
% comparaison des données de quelques profils aux bases de reference

if irun.VERIF_FLAG
    rep='VERIF_FLAG';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    %verif_flag_1profil(floatname,dacname,n_prof,'DEPTH_ZOOM',1000,'REFERENCE','ctd');
    verif_flag_1profil(floatname,dacname,n_prof,'DEPTH_ZOOM',1000, 'REFERENCE','argo');

    eval('cd ..');
    init_path('clear',rep,rdir);
end
%===========================%
% comparaison au profil de reference fait lors de la mise à l'eau
if irun.VERIF_PROF1_raw
    rep='VERIF_PROF1';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    VERIF_PROFIL1(floatname,dacname,1,config_campaign,'NB_PROF',10,'DIRECTION','A','ZOOM',1000);  % ici on compare le 1D avec le profil campagne, on trace les 10 premiers profils argo
    eval('cd ..');
    init_path('clear',rep,rdir);
end
%===========================%
% OPTIONNEL : comparaison des données Argo du flotteur a des profils (Argo) proches dans le temps et
% l'espace
if irun.FIND_CLOSE_float
    rep='COMPARE_FLOAT_REF_TR';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    FIND_close_floats(floatname,dacname,'UPDATE',1,'PRINT',1,'DATATYPE','raw')
    eval('cd ..');
    init_path('clear',rep,rdir);
end

%===========================%
% analyse OW
if irun.OW
    rep='LPO_CODES_ATLN_NEW';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    MAIN_dmqc_speed(floatname,dacname,numconfig_ow,'RECREATE_MAT',1);
    eval('cd ..');
    init_path('clear',rep,rdir);
    
    %===========================%
    % generation d'un rapport intermediaire & fichier config pour l'etape suivante (i.e. corrections)
    
    rep='DOC';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    generate_doc_overleaf(floatname,dacname,numconfig_ow,'SUBTITLE',['Float ' floatname ' (' name_campaign ')'],'PROFREF',n_prof)%
    eval('cd ..');
    init_path('clear',rep,rdir);
end

%===========================%
% Correction dans les fichiers
if irun.CORRECTIONS
    rep='CORRECTIONS';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    MAIN_write_dmqc_files(floatname,dacname) % A porter TR
    eval('cd ..');
    init_path('clear',rep,rdir);
end
%===========================%
%comparaison des données ajustées au profil de reference fait lors de la mise à l'eau
if irun.VERIF_PROF1_adj
    rep='VERIF_PROF1';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    VERIF_PROFIL1(floatname,dacname,1,config_campaign,'DIRECTION','A','NB_PROF',10,'DATATYPE','adj','DATAREP','DIR_DM_FILES','ZOOM',1000)
    eval('cd ..');
    init_path('clear',rep,rdir);
end
%===========================%
% on refait tourner les diag mais sur les données ajusté apres correction du CPCOR
if irun.PLOTDATA_adj
    rep='PROG_QC2015';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    %PLOTDATA_func(floatname,dacname,'DATATYPE','adj','DATAREP','DIR_DM_FILES','USEFLAG',1,'PRINT',1)   % attention matlab en display
    PLOTDATA_func(floatname,dacname,rdir,'DATATYPE','adj','DATAREP','DIR_DM_FILES','USEFLAG',1,'PRINT',1)   % attention matlab en display
    eval('cd ..');
    init_path('clear',rep,rdir);
end

%===========================%
% generation du rapport final
if irun.DOC
    rep='DOC';
    init_path('add',rep,rdir);
    eval(fullfile('cd ./',rep));
    %n_prof=0;% A revoir avec CC
    %generate_doc(floatname,dacname,numconfig_ow,'SUBTITLE',['Float ' floatname ' (bocats)'],'PROFREF',n_prof)%
    %generate_doc_overleaf(floatname,dacname,numconfig_ow,'SUBTITLE',['Float ' floatname ' (' name_campaign ')'],'PROFREF',n_prof,'PROFFLAG','1D')%
    generate_doc_overleaf(floatname,dacname,numconfig_ow,'SUBTITLE',['Float ' floatname ' (' name_campaign ')'],'PROFREF',n_prof)
    eval('cd ..');
    init_path('clear',rep,rdir);
    % uploaderles fichiers qui se trouvent dans ./DOC/OVERLEAF/ sur
    % https://www.overleaf.com/ et compiler le .tex
end
%===========================%
% sauvegarde de ce programme
if exist ('./TEMPLATES')==0
    mkdir ./TEMPLATES
end
copyfile ('MAIN_template.m', ['./TEMPLATES/MAIN_' floatname '.m'])



