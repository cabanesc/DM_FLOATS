%
% fct_pltmap
%
% function [hf,ha]=fct_pltmap(zone_visu,reso,proj)
%
% fonction qui trace une carte geographique d'une region
% donnee avec m_map
% la bathy est tracee et les contours de la terre
%
% INPUT :
% -------
% zone_visu = [latmin latmax lonmin lonmax]
% reso      = resolution de la bathy
%           = 'HR' haute resolution (AN_smithsandwell)
%           = 'LR' basse resolution (AN_smithsandwell_lowres)
% proj      = projection choisie
%           = 'lambert', 'mercator' ou autre (voir m_map)
%
%
%
% OUTPUT :
% --------
% hf,ha: handle pour la figure et les axes
%
%

function [hf,ha]=fct_pltmap(zone_visu,reso,proj)

if strcmp(reso,'HR') == 1
    if exist('D:\DATA\ARGO\data')
        load D:/DATA/ARGO/data/AN_smithsandwell
    else
        load /home5/pharos/argo/QCARGO/TRAITEMENT/data/AN_smithsandwell
    end
elseif strcmp(reso,'LR') == 1
  if exist('D:\DATA\ARGO\data')
        load D:/DATA/ARGO/data/AN_smithsandwell_lowres
    else
        load /home5/pharos/argo/QCARGO/TRAITEMENT/data/AN_smithsandwell_lowres 
  end
  bathy=bathy2;
  lon_bathy=lon_bathy2;
  lat_bathy=lat_bathy2;

  %load /home/revellata/vthierry/ATLNORD/ANE/HYDRO_ANE/ATLAS/ANE_smithsandwell.mat

 % load /home/revellata/vthierry/ATLNORD/BATHY/bathyhelen.mat
  %bathy=bathyhelen;
  %lon_bathy=atl3lon;
  %lat_bathy=atl3lat;
end

hf= figure ;
set(gca,'fontsize',10)
m_proj(proj,'long',[zone_visu(3) zone_visu(4)],...
                 'lat',[zone_visu(1) zone_visu(2)]);
m_grid('box','fancy','tickdir','in');	     

ha=gca;

hold on

v=[-2000 -2000];   
[n,p]=m_contour(lon_bathy,lat_bathy,bathy,v,'color',[.6 .6 .6]);
v=[-1000 -1000];
[q,r]=m_contour(lon_bathy,lat_bathy,bathy,v,'color',[.7 .7 .7]);
v=[-200 -200]; 
[c,h]=m_contour(lon_bathy,lat_bathy,bathy,v,'color',[.8 .8 .8]);
v=[0 0]; 
[c,h]=m_contour(lon_bathy,lat_bathy,bathy,v,'k','linewidth',2);

xlabel('Longitude')
ax=ylabel('Latitude');
getpo=get(ax,'position')
getpo(1) 
getpo(2)
getpo(1)=getpo(1) +getpo(1)*5/100;
set(ax,'position',getpo)
