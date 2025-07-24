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
% MAP_VISU = 
%           zone_visu = [latmin latmax lonmin lonmax]
%           reso      = resolution de la bathy
%                     = 'HR' haute resolution (AN_smithsandwell)
%                     = 'LR' basse resolution (AN_smithsandwell_lowres)
%           proj      = projection choisie
%                     = 'lambert', 'mercator' ou autre (voir m_map)
% PARAM =  
%          isarctic   =1 zone arctique
%                     =0 ailleurs
% OUTPUT :
% --------
% hf,ha: handle pour la figure et les axes
%
%

function [hf,ha]=fct_pltmap(MAP_VISU,PARAM)


hf= figure ;
set(hf,'Visible','on')
set(gca,'fontsize',18)


% Modified by T. Reynaud 08/09/2020

yt_min= MAP_VISU.min_lat;
yt_max= MAP_VISU.max_lat;
xt_min= MAP_VISU.min_lon;
xt_max= MAP_VISU.max_lon;

fact=0.1;
% m_proj(proj,'long',[floor(zone_visu(3)) ceil(zone_visu(4))],...
%                  'lat',[floor(zone_visu(1)) ceil(zone_visu(2))]);
%m_grid('box','fancy','tickdir','in');	     
%m_grid('tickdir','in');	

if PARAM.isarctic==0;
    m_proj( MAP_VISU.proj,'longitude',[xt_min xt_max],'latitude',[yt_min yt_max]);
    [bathy,lon_bathy,lat_bathy]=m_tbase([xt_min xt_max yt_min yt_max]);
    m_grid('color',[0 0 0],'linestyle',':');
    
else
    set(gca,'fontsize',12)
    m_proj(MAP_VISU.proj,'lat',(yt_max+yt_min)/2,'lon',(xt_min+xt_max)/2,'radius',min(yt_max-yt_min+1,50),'rectbox','off');
    [bathy,lon_bathy,lat_bathy]=m_tbase( [-180 180 0 90]);
    m_grid('color',[0 0 0],'linestyle',':');
    
end
             
%m_proj(proj,'longitude',[xt_min xt_max],'latitude',[yt_min yt_max]);
%m_grid('xtick',[xt_min:1/fact:xt_max],'ytick',[yt_min:1/fact:yt_max],'color',[0 0 0],'linestyle','-.');
m_grid('color',[0 0 0],'linestyle',':');


% Modified by T. Reynaud 07/09/2020
%m_proj(proj,'long',[floor(zone_visu(3)) ceil(zone_visu(4))],...
%                 'lat',[floor(zone_visu(1)) ceil(zone_visu(2))]);
%m_grid('box','fancy','tickdir','in');	     
%m_grid('tickdir','in');	     

ha=gca;

hold on

[bathy,lon_bathy,lat_bathy]=m_tbase([MAP_VISU.zone_visu(3) MAP_VISU.zone_visu(4) MAP_VISU.zone_visu(1) MAP_VISU.zone_visu(2)]);

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
