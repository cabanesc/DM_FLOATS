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

function [hf,ha]=fct_pltmap(zone_visu,reso,proj,park_press,prof_press)

screenSize = get(0, 'ScreenSize');

hf= figure('Position', [screenSize(3)*(1/4) screenSize(4)*(1/3) screenSize(3)*(1/2) screenSize(4)*(2/3)-90]) ;
%hf= figure('Position', [1 screenSize(4)*(1/3) screenSize(3) screenSize(4)*(2/3)-90]) ;
set(hf,'Visible','on')
set(gca,'fontsize',12)

% Modified by T. Reynaud 08/09/2020

fact=1;
yt_min= floor(zone_visu(1)*fact)/fact;
yt_max= ceil( zone_visu(2)*fact)/fact;
xt_min= floor(zone_visu(3)*fact)/fact;
xt_max= ceil( zone_visu(4)*fact)/fact;

% m_proj(proj,'long',[floor(zone_visu(3)) ceil(zone_visu(4))],...
%                  'lat',[floor(zone_visu(1)) ceil(zone_visu(2))]);
%m_grid('box','fancy','tickdir','in');	     
%m_grid('tickdir','in');	     
             
m_proj('miller','longitude',[xt_min xt_max],'latitude',[yt_min yt_max]);
           
ha=gca;

hold on

%[bathy,lon_bathy,lat_bathy]=m_tbase([zone_visu(3) zone_visu(4) zone_visu(1) zone_visu(2)]);
[bathy,lon_bathy,lat_bathy]=m_tbase([xt_min  xt_max yt_min yt_max]);

isoDerive = floor((park_press-30)/100)*100;
isoProfil1 = ceil((park_press+30)/100)*100;
isoProfil2 = floor((prof_press-30)/100)*100;
isoProfilLevels = [-isoProfil1:-200:-isoProfil2]';
if (length(isoProfilLevels) == 1)
      isoProfilLevels = [isoProfilLevels isoProfilLevels];
end


load ('bathy_28_colormap.mat');
newmap4=newmap3;
%newmap4=[newmap3;newmap3(end,:);newmap3(end,:);[1 0 1];[1 0 1];[1 0 1]];
cvec=[-12000:500:1000,100000];
[c,h]=m_contourf(lon_bathy,lat_bathy,bathy,cvec);
set(h,'LineStyle','None');
caxis([-10500 1000])
colormap(newmap4);
m_grid('xtick',[xt_min:1/fact:xt_max],'ytick',[yt_min:1/fact:yt_max],'color',[0 0 0],'linestyle',':');

% p = get(h,'Children');
% thechild=get(p,'CData');
% cdat=cell2mat(thechild);
% for i=1:length(cvec)-1
%     set(p(cdat>=cvec(i)& cdat< cvec(i+1)),'Facecolor',newmap4(i,:),'LineStyle','none')
% end

v=[-park_press-30 -park_press+30];
[q,r]=m_contour(lon_bathy,lat_bathy,bathy,v,'g');
v=[0:-200:-isoDerive];
[q,r]=m_contour(lon_bathy,lat_bathy,bathy,[0:-200:-isoDerive],'b');
clabel(q,r,v(1:4:end),'LabelSpacing',144*4,'FontSize',8,'FontWeight','normal','Color','b')

clabel(q,r,[1000 2000 3000])
v=[-prof_press-30 -prof_press+30];   
[n,p]=m_contour(lon_bathy,lat_bathy,bathy,v,'r');
[n,p]=m_contour(lon_bathy,lat_bathy,bathy,isoProfilLevels,'m');
clabel(n,p,isoProfilLevels(1:4:end),'LabelSpacing',144*4,'FontSize',8,'FontWeight','normal','Color','m')
v=[-200 -200]; 
%[c,h]=m_contour(lon_bathy,lat_bathy,bathy,v,'color',[.8 .8 .8]);
v=[0 0]; 
[c,h]=m_contour(lon_bathy,lat_bathy,bathy,v,'k','linewidth',2);

xlabel('Longitude')
ax=ylabel('Latitude');
getpo=get(ax,'position');

getpo(1)=getpo(1) +getpo(1)*5/100;
set(ax,'position',getpo)
