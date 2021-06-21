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


hf= figure ;
set(hf,'Visible','on')
set(gca,'fontsize',18)


% Modified by T. Reynaud 08/09/2020

fact=1;
yt_min= floor(zone_visu(1)*fact)/fact;
yt_max=  ceil( zone_visu(2)*fact)/fact;
xt_min= floor(zone_visu(3)*fact)/fact;
xt_max= ceil( zone_visu(4)*fact)/fact;

fact=0.1;
% m_proj(proj,'long',[floor(zone_visu(3)) ceil(zone_visu(4))],...
%                  'lat',[floor(zone_visu(1)) ceil(zone_visu(2))]);
%m_grid('box','fancy','tickdir','in');	     
%m_grid('tickdir','in');	     
             
m_proj(proj,'longitude',[xt_min xt_max],'latitude',[yt_min yt_max]);
m_grid('xtick',[xt_min:1/fact:xt_max],'ytick',[yt_min:1/fact:yt_max],'color',[0 0 0],'linestyle','-.');


% Modified by T. Reynaud 07/09/2020
%m_proj(proj,'long',[floor(zone_visu(3)) ceil(zone_visu(4))],...
%                 'lat',[floor(zone_visu(1)) ceil(zone_visu(2))]);
%m_grid('box','fancy','tickdir','in');	     
%m_grid('tickdir','in');	     

ha=gca;

hold on

[bathy,lon_bathy,lat_bathy]=m_tbase([zone_visu(3) zone_visu(4) zone_visu(1) zone_visu(2)]);

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
getpo=get(ax,'position');

getpo(1)=getpo(1) +getpo(1)*5/100;
set(ax,'position',getpo)
