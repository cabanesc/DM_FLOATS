# DM_FLOATS : How to install 

1. Clone the DM_FLOATS repo : git clone https://github.com/cabanesc/DM_FLOATS.git    => {REP_CODES}/DM_FLOATS

2. Add toolbox m_map1.4m in {REP_CODES}/DM_FLOATS/lib/  ; Change PATHNAME in m_etopo2.m l.47 (path to  ETOPO2v2g_i2_MSB.bin)

3. Dowload owc      https://github.com/ArgoDMQC/matlab_owc/releases/tag/v3.0.0   => {REP_CODES}/OWC/matlab_owc-3.0.0 

4. Download the last CTD reference database and the last ARGO reference database     => {REP_CODES}/OWC/matlab_ow-3.0.0/data/climatology/argo_profiles/
                                                                                     => {REP_CODES}/OWC/matlab_ow-3.0.0/data/climatology/historical_ctd/
5. run {REP_CODES}/DM_FLOATS/LPO_CODES_ATLN_NEW/genere_wmoboxes_matfile.m to generate {REP_CODES}/OWC/matlab_ow-3.0.0/data/constants/wmo_boxes_ctd.mat, wmo_boxes_argo.mat and wmo_boxes_ctdandargo.mat

6. Download a recent DOI of the full ARGO GDAC (needed to find nearby floats)
   http://www.argodatamgt.org/Access-to-data/Argo-DOI-Digital-Object-Identifier      => {REP_DATA}/ARGO_DOI/

7. Download the netcdf files of the float that is processed in Delayed Mode from the GDAC ftp server   =>{REP_DATA}/FTP_ORIG/
   ex: wget -r --user="anonymous" --password="your.email@institute.fr" ftp://ftp.ifremer.fr/ifremer/coriolis/argo/dac/coriolis/6901762/

8. Create a repository where you will store the netcdf files with intermediate corrections (FLAGS, CPCOR)   => {REP_DATA}/FTP_MODIF/

9. Create a repository where you will store the netcdf files with the final DM correction  (OWC adjustements) => {REP_DATA}/DM_FILES/

10. Create a repository where you will store plots from the various programs => {REP_DATA}/PLOTS/

11. Create a repository where you will store input/ouptut data from OWC (source, calib, mapped ,plot) => {REP_DATA}/data/

12. Download two topo files : topo.onetenthdeg.nc and topo.onedeg.nc   => {REP_DATA}/TOPO/

13. If a launch CTD is available, mat file can be put here  => {REP_DATA}/CAMPAIGN/ 

14. You will need pdflatex to create reports. As an alternative, it is possible to compile .tex online : overleaf  https://fr.overleaf.com/

15. Modify config_template.txt and save it to config.txt

16. Modify rdir in MAIN_template.m


