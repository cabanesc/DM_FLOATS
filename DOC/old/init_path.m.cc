% define_path.m
% Define path for the OW programs
% C.Cabanes 25/07/2012

% define path
close all

path(pathdef)

addpath('../lib/seawater_330_its90/')
addpath('../lib/libargo')
addpath('../lib/ext_lib')
addpath('../')
addpath('./util/')
addpath('/home1/homedir5/perso/ccabanes/dvlpRD/Argo/WorkOnBase/Trajectoire/Decodeur/soft/soft/m_map1.4e')
addpath('/home1/homedir5/perso/ccabanes/matlab/tool_Ocean/cmocean')
addpath('../PROG_QC2015/')
% % addpath('/home1/triagoz/matlab/outils_matlab/seawater/seawater_330_its90_lpo/')
% % addpath('/home1/triagoz/matlab/outils_matlab/lpo')
% % addpath('/home1/triagoz/matlab/outils_matlab/m_map1.4h')
% % addpath('/home1/homedir5/perso/ccabanes/dvlpRD/Argo/WorkOnBase/Trajectoire/Decodeur/soft/soft/m_map1.4e')
% % addpath('/home1/homedir5/perso/ccabanes/dvlpRD/Argo/Lib_forCathy/ObsInSitu/Lib_Argo/')
% % addpath('/home1/homedir5/perso/ccabanes/dvlpRD/Argo/Lib_forCathy/ObsInSitu/Lib_Argo/RWnetcdf/R2008b/')
% % addpath('/home1/homedir5/perso/ccabanes/dvlpRD/Argo/TD/CHECK_MOCCA/PROG_QC2015/')
% % REPERTOIRE des Donnees/Resultats  (./OW/data/)

% DIR_PLOT = ['/export/home1/ccabanes/data/QCARGO/TRAITEMENT/plot/'];
% DIR_DATA = ['/export/home1/ccabanes/data/QCARGO/TRAITEMENT/data/'];

% DIR_TRAJ = ['/home1/homedir5/perso/ccabanes/dvlpRD/Argo/WorkOnBase/Trajectoire/DISPLAY/'];

% %DIR_OW = ['/home1/homedir5/perso/ccabanes/dvlpRD/Argo/TD/CHECK_FLOAT/LPO_CODES_ATLN_NEW/data/float_plots/CONFIG'];
% DIR_PLOTOW = ['/export/home1/ccabanes/data/QCARGO/TRAITEMENT/data/float_plots/CONFIG'];

% DIR_DMQC=['/export/home1/ccabanes/data/DM_FILES/coriolis/'];  ==DIR_OUT => DIR_DM_FILES
% %creation des repertoires contenant les resultats (si ils n'existent pas)

% DIR_CODES=['/home1/homedir5/perso/ccabanes/dvlpRD/Argo/TD/CHECK_FLOATS/'];
% % REPERTOIRE GDAC FTP 

% %DIR_FTP=['/home/coriolis_exp/spool/co05/co0508/dac/'];
% DIR_FTP=['/export/home1/ccabanes/data/FTP_ARGO/'];



                    