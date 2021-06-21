function [hp,cb]=coliso(x,y,C,unit,lim,val,map);
%
% © B.Ferron (11-96), bferron@ifremer.fr
%
% Name:
% COLISO is a contour plot that fills each contour
% 	level with a solid color.
%
% Synopsis and description:
% COLISO(x,y,C,unit,iso) 
%	x,y are vectors and C is the corresponding matrix
%	to be contoured: size(C)=[length(y) length(x)]
%	unit is a string expression (styled text) which is displayed
%	on the colorbar to specify the unit of the values of C.
%	Type help stext to see way styled text objects are created
%	(e.g if C is a matrix of temperature in Celcius then 
%	unit='{\degrees}C' displays °C on the colorbar)
%	iso is the vector wich contains the contours
%	(e.g iso=-4000:500:-3000 means that you want to contour the
%	values of C from -4000 to -3000 changing the solid color
%	every 500)
% COLISO(x,y,C,unit,iso,map) 
%	use the RGB colormap map instead of the default colormap(jet)
% COLISO(x,y,C,unit,iso,map,'horiz') 
%	same command with an horizontal color scale
% COLISO(x,y,C,unit,iso,val,map)
%	Allow you to plot the value val of the matrix C with the
%	RGB color map. The colors between isocontours are set to 
%	the colormap(jet). This function may be useful if you have
%	a matrix C with NaN values and masked values.
% [HP,HC]=COLISO(...) return the handle to the pcolor HP and the colorbar HC
%


axis([min(x) max(x) min(y) max(y)])

cmin=min(lim);
cmax=max(lim);
diso=lim(2)-lim(1);
iso=cmin:diso:cmax;
niso=length(iso)-1;
[nl,nc]=size(C);
ind=1:nl*nc;

if nargin==7 & ~ischar(map)
   tmp0=find(C==val);
   ind(tmp0)=[];
end

tmp=find(C(ind)<cmin);
C(ind(tmp))=cmin*ones(size(tmp));
tmp=find(C(ind)>cmax);
C(ind(tmp))=(cmax-diso/2)*ones(size(tmp));

colormap(jet(niso))
map0=colormap;


if nargin==7 & ~ischar(map)

   C(tmp0)=(cmax+diso/2)*ones(size(tmp0));
   z=fix(1+(C-cmin)/diso);
   caxis([1 niso+1])
   colormap([map0;map])
   cb=colorbar;
   pos=get(gca,'position');
   rat=get(gca,'aspect');
   poscb=get(cb,'pos');
%   if isnan(rat)~=[1 1]
%      set(cb,'position',[poscb(1) pos(2)+pos(4)/2-pos(3)/(2*rat(1)) poscb(3) pos(3)/rat(1)])
%   end
   
   hold on
   hp=pcolor(x,y,z');
   set(hp,'edgecolor','none');
   set(hp,'zdata',-ones(size(get(hp,'cdata'))))

   chil=get(cb,'child');
   nch=length(chil);
   set(chil(nch),'cdata',[1:niso]')
   cby=get(cb,'ylim');
   posc=get(chil(1),'posi');
   for i=1:length(chil)-1
       set(chil(i),'string','')
   end  
   dy=(cby(2)-cby(1))/niso;
   ca=gca;
   axes(cb)
   pause(1)
   for i=1:niso+1
       text(posc(1)-0.3,cby(1)+(i-1)*dy,[' ' num2str(iso(i)) ' ' unit])
   end  
   set(cb,'ytick',[],'xtick',[])
   axes(ca)
else

   if nargin==6 | ischar(map)
      colormap(val)
      map0=colormap;
   end
   z=fix(1+(C-cmin)/diso);
   caxis([1 niso+1])
   colormap(map0)
   if exist('map')
      cb=colorbar('horiz');
   else
      cb=colorbar;
   end
% 
%   pos=get(gca,'position');
%   rat=get(gca,'aspect');
%   poscb=get(cb,'pos');
%    
%   if isnan(rat)~=[1 1]
%      set(cb,'position',[poscb(1) pos(2)+pos(4)/2-pos(3)/(2*rat(1)) poscb(3) pos(3)/rat(1)])
%   end

   hold on

   
% Version 4.2c
%   chil=get(cb,'child');
%   nch=length(chil);
%   set(chil(nch),'cdata',[1:niso]')
%   cby=get(cb,'ylim');
%   posc=get(chil(1),'posi');
%   for i=1:length(chil)-1
%       set(chil(i),'string','')
%   end  
%   dy=(cby(2)-cby(1))/niso;
   ca=gca;
%   axes(cb);
%   pause(1)
%   for i=1:niso+1
%       text(posc(1)-0.3,cby(1)+(i-1)*dy,...
%	    [num2str(iso(i)) ' ' unit])
%	text(posc(1)-0.3,cby(1)+(i-1)*dy,[ num2str(iso(i)) ' ' unit])
%   end
   if exist('map')
      xt=get(cb,'xticklabel');
      set(cb,'xticklabel',[int2str(iso(str2num(xt))') char(ones(length(xt),1)*unit)])
   else
      yt=get(cb,'yticklabel');
      set(cb,'yticklabel',[int2str(iso(str2num(yt))') char(ones(length(yt),1)*unit)])
%   set(cb,'ytick',[],'xtick',[])
   end
   axes(ca)

   hp=pcolor(x,y,z');
   set(hp,'edgecolor','none');
   set(hp,'zdata',-ones(size(get(hp,'cdata'))))
end  
