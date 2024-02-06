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

function fct_pltmap_new(zone_visu,reso,proj)

% if strcmp(reso,'HR') == 1
%     if exist('D:\DATA\ARGO\data')
%         load D:/DATA/ARGO/data/AN_smithsandwell
%     else
%         load /home1/revellata/vthierry/BATHYMETRY/AN_smithsandwell
%     end
% elseif strcmp(reso,'LR') == 1
%   if exist('D:\DATA\ARGO\data')
%         load D:/DATA/ARGO/data/AN_smithsandwell_lowres
%     else
%         load /home1/revellata/vthierry/BATHYMETRY/AN_smithsandwell_lowres 
%   end
%   bathy=bathy2;
%   lon_bathy=lon_bathy2;
%   lat_bathy=lat_bathy2;
% 
%   %load /home/revellata/vthierry/ATLNORD/ANE/HYDRO_ANE/ATLAS/ANE_smithsandwell.mat
% 
%  % load /home/revellata/vthierry/ATLNORD/BATHY/bathyhelen.mat
%   %bathy=bathyhelen;
%   %lon_bathy=atl3lon;
%   %lat_bathy=atl3lat;
% end


%set(gca,'fontsize',12)
m_proj(proj,'long',[zone_visu(3)-1 zone_visu(4)+1],...
                 'lat',[zone_visu(1)-1 zone_visu(2)+1]);
warning off
m_grid('box','fancy','tickdir','out', 'linestyle', 'none','xtick',3);	     
warning on
ha=gca;

hold on

%[bathy,lon_bathy,lat_bathy]=m_tbase([zone_visu(3) zone_visu(4) zone_visu(1) zone_visu(2)]);
%keyboard
%[contourMatrix, contourHdl] = m_etopo2('contour', [-4000 -4000], 'color',[0 0 .5]);
%[contourMatrix, contourHdl] = m_etopo2('contour', [-3000 -3000], 'color',[.4 .4 .9]);
[contourMatrix, contourHdl] = m_etopo2('contour', [-2000 -2000], 'color',[.7 .7 .9]);
[contourMatrix, contourHdl] = m_etopo2('contour', [-1000 -1000], 'color',[.4 .9 .9]);
[contourMatrix, contourHdl] = m_etopo2('contour', [-200 -200], 'color',[.7 .9 .9]);
[contourMatrix, contourHdl] = m_etopo2('contour', [0 0], 'color',[0 0 0],'linewidth',2);
%  v=[-4000 -4000];   
%  [n,p]=m_contour(lon_bathy,lat_bathy,bathy,v,'color',[.3 .3 .3]);
%  v=[-3000 -3000];   
%  [n,p]=m_contour(lon_bathy,lat_bathy,bathy,v,'color',[.5 .5 .7]);
%  v=[-2000 -2000];   
%  [n,p]=m_contour(lon_bathy,lat_bathy,bathy,v,'color',[.6 .5 .6]);
%  v=[-1000 -1000];
%  [q,r]=m_contour(lon_bathy,lat_bathy,bathy,v,'color',[.7 .6 .7]);
%  v=[-200 -200]; 
%  [c,h]=m_contour(lon_bathy,lat_bathy,bathy,v,'color',[.8 .8 .8]);
%  v=[0 0]; 
%  [c,h]=m_contour(lon_bathy,lat_bathy,bathy,v,'k','linewidth',2);

xlabel('Longitude')
ax=ylabel('Latitude');
getpo=get(ax,'position');


getpo(1)=getpo(1) +getpo(1)*5/100;
set(ax,'position',getpo)
