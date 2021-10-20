%   General configuration file

% GENERAL
% -----------------

% miror of the GDAC FTP server (original floats profiles to process in DM, won't be modified)
DIR_FTP_CORIOLIS=/home/coriolis_exp/spool/co05/co0508/dac/

% Last DOI of the full ARGO GDAC (needed to FIND_CLOSE_FLOATS)
DIR_REFDATA_CORIOLIS=/home5/pharos/REFERENCE_DATA/ARGO/DATA/202106-ArgoData/dac/


% local directory to store float netcdf files you are working on (copied from DIR_FTP_CORIOLIS), will be modified (RT flags)
DIR_FTP=/export/home1/ccabanes/data/FTP_ARGO/

% Final files with DMQC corrections are put in this directory
DIR_DM_FILES=/export/home1/ccabanes/data/DM_FILES/


% place to save your plots from various codes 
DIR_PLOT=/export/home1/ccabanes/data/QCARGO/TRAITEMENT/plot/

% place for your input/ouptut data from OWC (source, calib, mapped ,plot)
DIR_DATA=/export/home1/ccabanes/data/QCARGO/TRAITEMENT/data/

% full path of the DM_FLOATS package directory
DIR_CODES=/home1/homedir5/perso/ccabanes/dvlpRD/Argo/TD/git/DM_FLOATS/

% OW codes and reference data are stored on a common disk (pharos)
%-----------------------------------------------------------------
% OWC path on the common disk (pharos)
DIR_OWC=/home5/pharos/argo/DMARGO/OW/

% OWC directory: Full path is [DIR_OWC VERSION_OWC] e.g /home5/pharos/argo/DMARGO/OW/matlab_ow-3.0.0_a/
VERSION_OWC=matlab_ow-3.0.0_a

% Reference data are in /home5/pharos/argo/DMARGO/OW/matlab_ow-3.0.0_a/data/climatology/
% run ./LPO_CODES_ATLN_NEW/genere_wmoboxes_matfile.m to generate /home5/pharos/argo/DMARGO/OW/matlab_ow-3.0.0_a/data/constants/wmo_boxes_ctd.mat, wmo_boxes_argo.mat and wmo_boxes_ctdandargo.mat

% Bathymetry files:
FILE_TOPO=/home5/pharos/argo/DMARGO/OW/TOPO/topo.onetenthdeg.nc
FILE_TOPO_LOW=/home5/pharos/argo/DMARGO/OW/TOPO/topo.onedeg.nc

% Reference data files. Full path is [DIR_OW DIR_OWC_CTD] e.g. /home5/pharos/argo/DMARGO/OW/LIMATOLOGY/CTD_for_DMQC_2019V01/
%DIR_OWC_CTD=CLIMATOLOGY/CTD_for_DMQC_2021V01/
%DIR_OWC_ARGO=CLIMATOLOGY/ARGO_for_DMQC_2020V01/
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
