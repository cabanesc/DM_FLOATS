%==========================================================================
% Auteur: DAVID Nicolas - IFREMER/ENSIETA (Adapted from create_float_source.mat by CC.)
% Date: 09/08/2007
% Objectif: Creer les fichiers sources .mat des campagnes oceanographiques
% a� partir des fichiers .netcdf.
% Modification:
%==========================================================================

%addpath('../matlab_codes/')

%po_system_configuration = load_configuration( ['config_ovide08_thalassa.txt'] );
%po_system_configuration = load_configuration( 'config_rrex17.txt');
%po_system_configuration = load_configuration( 'config_geovide.txt');
%po_system_configuration = load_configuration( 'config_bocats.txt');
po_system_configuration = load_configuration( 'config_catarina.txt');

%--------------------------------------------------------------------------
% Recherche des donnees
%--------------------------------------------------------------------------

netcdf_filename=[po_system_configuration.CAMPAGNE_DIRECTORY po_system_configuration.CAMPAGNE_FILE]
nc=netcdf.open(netcdf_filename,'nowrite');

latb     = netcdf.getVar(nc,netcdf.inqVarID(nc,'LATITUDE_BEGIN'));
lonb     = netcdf.getVar(nc,netcdf.inqVarID(nc,'LONGITUDE_BEGIN'));
juldb    = netcdf.getVar(nc,netcdf.inqVarID(nc,'JULD_BEGIN'));
latb(latb  ==-9999)   = NaN;
lonb(lonb  ==-9999)   = NaN;
juldb(juldb==-9999)   = NaN;

late     = netcdf.getVar(nc,netcdf.inqVarID(nc,'LATITUDE_END'));
lone     = netcdf.getVar(nc,netcdf.inqVarID(nc,'LONGITUDE_END'));
julde    = netcdf.getVar(nc,netcdf.inqVarID(nc,'JULD_END'));
late(late==-9999)     = NaN;
lone(lone==-9999)     = NaN;
julde(julde==-9999)   = NaN;

LAT   = mynanmean([latb late]')';
LONG  = mynanmean([lonb lone]')';
DATES = mynanmean([juldb julde]')';

%tabvar=['PRES';'TEMP';'PSAL'];

PRES     = netcdf.getVar(nc,netcdf.inqVarID(nc,'PRES'))';
fillvalue=netcdf.getAtt(nc,(netcdf.inqVarID(nc,'PRES')),'_FillValue');
ipres=find(PRES==fillvalue);
PRES(ipres) = NaN;

TEMP     = netcdf.getVar(nc,netcdf.inqVarID(nc,'TEMP'))';
fillvalue=netcdf.getAtt(nc,(netcdf.inqVarID(nc,'TEMP')),'_FillValue');
itemp=find(TEMP==fillvalue);
TEMP(itemp) = NaN;

SAL     = netcdf.getVar(nc,netcdf.inqVarID(nc,'PSAL'))';
fillvalue=netcdf.getAtt(nc,(netcdf.inqVarID(nc,'PSAL')),'_FillValue');
ipsal=find(SAL==fillvalue);
SAL(ipsal) = NaN;



%tabvar=['PRES';'TEMP';'PSAL'];
%for ivar=1:3
%	tempo=netcdf.getVar(nc,netcdf.inqVarID(nc,tabvar(ivar,:)(:)));
%	fillvalue=nc{tabvar(ivar,:)}.FillValue_(:);
%	tempo(find(tempo==fillvalue))=NaN;
%	eval([lower(tabvar(ivar,:)) '=tempo;']);
%end

netcdf.close(nc);
    
clear netcdf_filename nc latb lonb juldb late lone julde fillvalue 
    

%--------------------------------------------------------------------------
% Ecriture du fichier .mat
%--------------------------------------------------------------------------

mat_filename=[po_system_configuration.DATA_DIRECTORY strcat(po_system_configuration.CAMPAGNE_MAT, po_system_configuration.POSTFIX ) ]

save(mat_filename,'DATES','LAT','LONG','PRES','TEMP','SAL')

clear all

