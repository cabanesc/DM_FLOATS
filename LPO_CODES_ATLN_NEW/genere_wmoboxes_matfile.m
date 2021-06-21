% Genere les fichiers wmo_boxes_xxx.mat en fonction
% de l'existence des fichiers boxnum.mat

inipath='/home/lpoargo1/DMARGO/OW/VERSION_1_1_3/data/';

% Chargement du fichier initial, r�cup�ration du tableau la_wmo_boxes
load([inipath 'constants/wmo_boxes.mat']);
la_wmo_boxes_ini=la_wmo_boxes;

% Fichier contenant les donn�es CTD seulement pour les bo�tes WMO
% existantes
la_wmo_boxes(:,2)=0;
la_wmo_boxes(:,4)=0;
[nwmo,ncol]=size(la_wmo_boxes);
for ind=1:nwmo
    if exist([inipath 'climatology/historical_ctd/ctd_' num2str(la_wmo_boxes(ind,1)) '.mat']) == 2
        la_wmo_boxes(ind,2)=1;
    end
end
la_wmo_boxes_ctd=la_wmo_boxes(:,2);
eval(['save ' inipath 'constants/wmo_boxes_ctd.mat  la_wmo_boxes']);

% Fichier contenant les donn�es CTD seulement de toutes les bo�tes WMO 
la_wmo_boxes(:,2)=1;
eval(['save ' inipath 'constants/wmo_boxes_allctd.mat  la_wmo_boxes']); 

% Fichier contenant les donn�es ARGO seulement de toutes les bo�tes WMO 
la_wmo_boxes(:,2)=0;
la_wmo_boxes(:,4)=1;
eval(['save ' inipath 'constants/wmo_boxes_allargo.mat  la_wmo_boxes']); 


% Fichier contenant les donn�es ARGO seulement pour les bo�tes WMO
% existantes
la_wmo_boxes(:,2)=0;
la_wmo_boxes(:,4)=0;
[nwmo,ncol]=size(la_wmo_boxes);
for ind=1:nwmo
    if exist([inipath 'climatology/argo_profiles/argo_' num2str(la_wmo_boxes(ind,1)) '.mat']) == 2
        la_wmo_boxes(ind,4)=1;
    end
end
la_wmo_boxes_argo=la_wmo_boxes(:,4);
eval(['save ' inipath 'constants/wmo_boxes_argo.mat  la_wmo_boxes']); 
   
% Fichier contenant les donn�es ARGO et les donn�es CTD pour les bo�tes WMO
% existantes
la_wmo_boxes(:,2)=la_wmo_boxes_ctd;
la_wmo_boxes(:,4)=la_wmo_boxes_argo;
eval(['save ' inipath 'constants/wmo_boxes_ctdandargo.mat  la_wmo_boxes']); 

% AJOUT permettant d'acclerer la lecture des données de reference : fichier ne contenant que lat et
% lon pour chaque boite
if exist([inipath 'climatology/historical_ctd/lonlat'])==0
    mkdir([inipath 'climatology/historical_ctd/lonlat'])
end
% genere les fichiers lat long date pour chaque boite.
for ind=1:nwmo
    file=[inipath 'climatology/historical_ctd/ctd_' num2str(la_wmo_boxes(ind,1)) '.mat'];
    fileout=[inipath 'climatology/historical_ctd/lonlat/ctd_' num2str(la_wmo_boxes(ind,1)) '.mat'];
    if exist([inipath 'climatology/historical_ctd/ctd_' num2str(la_wmo_boxes(ind,1)) '.mat']) == 2
        l=load(file);
        l=rmfield(l,{'pres','sal','temp','ptmp'});
        lat=l.lat;
        long=l.long;
        dates=l.dates;
        source=l.source;
        
        save(fileout,'lat','long','dates','source')
    end
end

if exist([inipath 'climatology/argo_profiles/lonlat'])==0
    mkdir([inipath 'climatology/argo_profiles/lonlat'])
end

for ind=1:nwmo
    file=[inipath 'climatology/argo_profiles/argo_' num2str(la_wmo_boxes(ind,1)) '.mat'];
    fileout=[inipath 'climatology/argo_profiles/lonlat/argo_' num2str(la_wmo_boxes(ind,1)) '.mat'];
    if exist([inipath 'climatology/argo_profiles/argo_' num2str(la_wmo_boxes(ind,1)) '.mat']) == 2
        l=load(file);
        l=rmfield(l,{'pres','sal','temp','ptmp'});
        lat=l.lat;
        long=l.long;
        dates=l.dates;
        source=l.source;
        
        save(fileout,'lat','long','dates','source')
    end
end
