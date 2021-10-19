function [gl_lat,gl_lon,gl_cc,gl_day,gl_mon,gl_year,gl_sal,gl_temp,gl_oxy,gl_z] = extract_hydro(latdeb, latfin, londeb, lonfin)
% Cette fonction permet l'extraction de stations hydrologiques
% de la base de donnees de Thierry Reynaud
%
% Les parametres d'entrees sont:
% latdeb= borne sud du domaine
% latfin= borne nord du domaine
% londeb= borne ouest du domaine
% lonfin= borne est du domaine
%
% Les parametres de sorties sont:
% gl_lat: latitudes des stations
% gl_lon: longitudes des stations
% gl_cc:  Code pays NODC des stations
% gl_day/gl_mon/gl_year:  Date (jour/mois/annee) des stations
% gl_z:    immersions des stations
% gl_sal:  salinites [psu] des stations
% gl_temp: temperatures insitu [Deg. Celcius]des stations
% gl_oxy:  oxygene [umol/kg] des stations
%

path(path,'/home/mengam/treynaud/MATLAB/TOOLS/');

%Creation des matrices de stockage
gl_lat=[];
gl_lon=[];
gl_day=[];
gl_mon=[];
gl_year=[];
gl_sal=[];
gl_temp=[];
gl_oxy=[];
gl_z=[];
gl_cc=[];

file='/home/goodhope/treynaud/ANALYZED_DATA/STANDARD/FILTERED/filtered_good_final.v2005.bin';

% Ouverture du fichier BIMG: ieee-sun c'est du ieee big-endian
fid = fopen(file, 'r','ieee-be');

%Lecture du nombre de stations hydrologiques contenues dans la base:
%Lecture du nombre de niveaux verticaux:
nbid = fread(fid, 1, 'int32');
tab = fread(fid, 2, 'int32');
nbid = fread(fid, 1, 'int32');

nmax=tab(1);nz=tab(2);

nbid = fread(fid, 1, 'int32');
tab = fread(fid, nz, 'float32');
nbid = fread(fid, 1, 'int32');
gl_z=tab;

clear tab;
nbid = fread(fid, 1, 'int32');
tab = fread(fid, [2 nmax], 'char');
nbid = fread(fid, 1, 'int32');
gl_cc = char(tab)';

clear tab;
nbid = fread(fid, 1, 'int32');
tab = fread(fid, nmax, 'int32');
nbid = fread(fid, 1, 'int32');
gl_day=tab;

clear tab;
nbid = fread(fid, 1, 'int32');
tab = fread(fid, nmax, 'int32');
nbid = fread(fid, 1, 'int32');
gl_mon=tab;

clear tab;
nbid = fread(fid, 1, 'int32');
tab = fread(fid, nmax, 'int32');
nbid = fread(fid, 1, 'int32');
gl_year=tab;

clear tab;
nbid = fread(fid, 1, 'int32');
tab = fread(fid, nmax, 'float32');
nbid = fread(fid, 1, 'int32');
gl_lat=tab;

clear tab;
nbid = fread(fid, 1, 'int32');
tab = fread(fid, nmax, 'float32');
nbid = fread(fid, 1, 'int32');
gl_lon=tab;

idx=find( gl_lat <= latfin & gl_lat >= latdeb &... 
gl_lon <= lonfin & gl_lon >= londeb);

gl_lat=gl_lat(idx);
gl_lon=gl_lon(idx);
gl_day=gl_day(idx);
gl_mon=gl_mon(idx);
gl_year=gl_year(idx);
gl_cc=gl_cc(idx,:);


clear tab;
nbid = fread(fid, 1, 'int32');
tab = fread(fid, [nz,nmax], 'float32');
nbid = fread(fid, 1, 'int32');
gl_temp=tab(:,idx);clear tab;


nbid = fread(fid, 1, 'int32');
tab = fread(fid, [nz,nmax], 'float32');
nbid = fread(fid, 1, 'int32');
gl_sal=tab(:,idx);clear tab;


nbid = fread(fid, 1, 'int32');
tab = fread(fid, [nz,nmax], 'float32');
nbid = fread(fid, 1, 'int32');
gl_oxy=tab(:,idx);clear tab;
