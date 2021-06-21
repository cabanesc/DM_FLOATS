%
%   ***************
%   *   Palette   *
%   ***************
%
% Appel : palette(cmin,cmax,nbre_couleurs,flag,flag2)
%
% Fabrique une palette de couleurs "propre", ou les limites 
% des couleurs correspondent aux labels.
%
%   -> cmin est la valeur minimum du champ a tracer
%   -> cmax               maximum
%   -> nbre_couleurs est le nombre de couleurs que l'on
%      souhaite dans la palette.
%   -> flag=1 si trace de la colorbar, 0 sinon
%   -> flag2 est la taille de la fonte de la colorbar 12, 14, 18...
%
% © Laurent DUBUS (11-1997), merci Bruno !
%

function [h]=palette3(cmin,cmax,nconts,flag,flag2)

caxis([cmin cmax]);
colormap(jet(nconts));
map=colormap;
set(gca,'colororder',map)
if (flag==1)
  h=colorbar('horiz');
  set(h,'YTick',(cmin:(cmax-cmin)/nconts*2:cmax),'FontSize',flag2)
end
