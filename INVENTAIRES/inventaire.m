rdir='/home1/homedir5/perso/ccabanes/dvlpRD/Argo/TD/DM_FLOATS_TR/'; % cc
addpath(fullfile(rdir,'/lib/libargo'));
addpath(fullfile(rdir,'/lib/ext_lib'));

CONFIG=load_configuration('config.txt');

index_file_name=[CONFIG.DIR_FTP_CORIOLIS '/../ar_index_global_prof.txt'];
index_data = read_index_prof(index_file_name);
    
[liste_profile_files,liste_floats_files] = select_floats_from_index(index_data,'centre',{'IF'});
cvs_pathfilename=['inventaire_co_' date '.csv'];
warning off
make_index_DM(CONFIG.DIR_FTP_CORIOLIS,liste_floats_files,cvs_pathfilename)

% MOCCA floats (BODC&ifremer)
liste_floats_files={
'coriolis/3901871/3901871_prof.nc'
'coriolis/3901901/3901901_prof.nc'
'coriolis/3901902/3901902_prof.nc'
'coriolis/3901903/3901903_prof.nc'
'coriolis/3901918/3901918_prof.nc'
'coriolis/3901919/3901919_prof.nc'
'coriolis/3901920/3901920_prof.nc'
'coriolis/3901921/3901921_prof.nc'
'coriolis/3901922/3901922_prof.nc'
'coriolis/3901923/3901923_prof.nc'
'coriolis/3901924/3901924_prof.nc'
'coriolis/3901925/3901925_prof.nc'
'coriolis/3901926/3901926_prof.nc'
'coriolis/3901927/3901927_prof.nc'
'coriolis/3901928/3901928_prof.nc'
'coriolis/3901929/3901929_prof.nc'
'coriolis/3901930/3901930_prof.nc'
'coriolis/3901931/3901931_prof.nc'
'coriolis/3901932/3901932_prof.nc'
'coriolis/3901933/3901933_prof.nc'
'coriolis/3901934/3901934_prof.nc'
'coriolis/3901935/3901935_prof.nc'
'coriolis/3901936/3901936_prof.nc'
'coriolis/3901937/3901937_prof.nc'
'bodc/3901943/3901943_prof.nc'
'bodc/3901944/3901944_prof.nc'
'bodc/3901945/3901945_prof.nc'
'bodc/3901951/3901951_prof.nc'
'bodc/3901954/3901954_prof.nc'
'bodc/3901955/3901955_prof.nc'
'bodc/3901956/3901956_prof.nc'
'bodc/3901964/3901964_prof.nc'
'bodc/3901965/3901965_prof.nc'
'bodc/3901970/3901970_prof.nc'
'bodc/3901971/3901971_prof.nc'
'bodc/3901972/3901972_prof.nc'
'bodc/3901980/3901980_prof.nc'
'bodc/3901982/3901982_prof.nc'
'bodc/3901983/3901983_prof.nc'
'bodc/3901984/3901984_prof.nc'}
cvs_pathfilename=['inventaire_MOCCA_' date '.csv'];
make_index_DM(CONFIG.DIR_FTP_CORIOLIS,liste_floats_files,cvs_pathfilename)
warning on