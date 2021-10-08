%   General configuration file

% GENERAL
% -----------------

% miror of the GDAC FTP server (original floats profiles, won't be modified)
DIR_FTP_CORIOLIS=/home/coriolis_exp/spool/co05/co0508/dac/

% Last DOI of the full ARGO GDAC
DIR_REFDATA_CORIOLIS=/home5/pharos/REFERENCE_DATA/ARGO/DATA/202106-ArgoData/dac/


% local directory to store float netcdf files you are working on (copied from DIR_FTP_CORIOLIS)
DIR_FTP=/export/home1/ccabanes/data/FTP_ARGO/

% output files with DMQC corrections are put in this directory
DIR_DM_FILES=/export/home1/ccabanes/data/DM_FILES/


% place to save your plots from various codes 
DIR_PLOT=/export/home1/ccabanes/data/QCARGO/TRAITEMENT/plot/

% place for your input/ouptut data from OWC (source, calib, mapped ,plot)
DIR_DATA=/export/home1/ccabanes/data/QCARGO/TRAITEMENT/data/

% full path of the DM_FLOATS directory
DIR_CODES=/home1/homedir5/perso/ccabanes/dvlpRD/Argo/TD/git/DM_FLOATS/

% OW codes and reference data are stored on a common disk (pharos)
%-----------------------------------------------------------------
% OWC path on the common disk (pharos)
DIR_OWC=/home5/pharos/argo/DMARGO/OW/

% OWC version to use  ('VERSION_LAST' or 'VERSION_1_1_9' or 'VERSION_2_1_0' ... see DIR_OWC for available versions)
VERSION_OWC=VERSION_3_0_1

% Bathymetry files:
FILE_TOPO=/home5/pharos/argo/DMARGO/OW/TOPO/topo.onetenthdeg.nc
FILE_TOPO_LOW=/home5/pharos/argo/DMARGO/OW/TOPO/topo.onedeg.nc

% Reference data files
DIR_OWC_CTD=CLIMATOLOGY/CTD_for_DMQC_2019V01/
DIR_OWC_ARGO=CLIMATOLOGY/ARGO_for_DMQC_2019V03/
% /Volumes/argo/DMARGO/OW/VERSION_2_1_0/data/climatology/historical_ctd/ctd_7503 ... pas en smb
% /Volumes/argo/DMARGO/OW/CLIMATOLOGY/CTD_for_DMQC_2019V01/ctd_7503


% CTD MADE AT float launch
% -------------
% CAMPAIGN CTD data are stored on a common disk (pharos)
DIR_CAMPAIGN=/home5/pharos/argo/QCARGO/TRAITEMENT/data/

% DM OPERATOR INFO
%-----------------
OPERATOR_ORCID_ID = https://orcid.org/0000-0003-3873-5259
OPERATOR_NAME = C. Cabanes
OPERATOR_INSTITUTION = LOPS/Ifremer
