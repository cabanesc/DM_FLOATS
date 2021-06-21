function [thetitle]=plot_bathy_glob_nofill(region,Topo)
% -========================================================
%   USAGE : plot_traj(FL)
%   PURPOSE : plot the trajectory of the float, on a map with bathy
% -----------------------------------
%   INPUT :
%     FL   (structure)  Float data
%   optional input :
%     Vec 1*n_cycle parameter to plot along the trajectory
%   OUTPUT :
%    thetitle (char) title of the plot
% -----------------------------------
%   HISTORY  : created (2009) ccabanes
%            : modified (yyyy) byxxx
%   CALLED SUBROUTINES: none
% ========================================================



if isfield(Topo,'lat')==0&isfield(Topo,'grd_lat')==1
   Topo.lat=Topo.grd_lat;
end
if isfield(Topo,'lon')==0&isfield(Topo,'grd_lon')==1
   Topo.lon=Topo.grd_lon;
end

Topo.lon.data=Topo.lon.data(1:end-1);
Topo.topo.data=Topo.topo.data(:,1:end-1);
Topo=shiftEW(Topo,'lon','grwch');
ii=find(Topo.lat.data>=region.latmin-1&Topo.lat.data<=region.latmax+1);
jj=find(Topo.lon.data>=region.lonmin-1&Topo.lon.data<=region.lonmax+1);
Topo.lon.data=Topo.lon.data(jj);
Topo.lat.data=Topo.lat.data(ii);
Topo.topo.data=Topo.topo.data(ii,jj);
[junk,isort]=sort(Topo.lon.data);
Topo.lon.data=Topo.lon.data(isort);
Topo.topo.data=Topo.topo.data(:,isort);

% load colormap
load ('colormap_bathy6.mat');


        
%keyboard
if min(size(Topo.lat.data))==1
    Topo.lat.data =repmat(Topo.lat.data,[1,length(Topo.lon.data)]);
    Topo.lon.data =repmat(Topo.lon.data',[size(Topo.lat.data,1),1]);
end

iy = (Topo.lat.data(:,1)>= region.latmin & Topo.lat.data(:,1)<= region.latmax);
ix = (Topo.lon.data(1,:)>= region.lonmin & Topo.lon.data(1,:)<= region.lonmax);
%  r=reshape(Topo.lat.data(iy,ix),sum(ix)*sum(iy),1);
%  t=reshape(Topo.lon.data(iy,ix),sum(ix)*sum(iy),1);
%  y=reshape(Topo.topo.data(iy,ix),sum(ix)*sum(iy),1);
%  y(y>0)=NaN;



box on
grid on
hold on
% if isempty(regionI)||(region.lonmin < regionI.lonmin|region.lonmax>regionI.lonmax|region.latmin<regionI.latmin|region.latmax>regionI.latmax)

cvec=[-4000,-2000,-1000,0,1000,100000];

for i=1:length(cvec)-1
[c,h]=contour(Topo.lon.data(iy,ix),Topo.lat.data(iy,ix),Topo.topo.data(iy,ix),[cvec(i) cvec(i)],'Color',newmap(i,:));
end
%[c,h]=contourf(Topo.lon.data(iy,ix),Topo.lat.data(iy,ix),Topo.topo.data(iy,ix),[cvec(15:end)]);

[c,h]=contourf(Topo.lon.data(iy,ix),Topo.lat.data(iy,ix),Topo.topo.data(iy,ix),[cvec([4:end])]);
%keyboard
p = get(h,'Children');
thechild=get(p,'CData');
%keyboard
if length(thechild)>1
cdat=cell2mat(thechild);
else
cdat=thechild;
end
for i=[4:length(cvec)-1]
    set(p(cdat>=cvec(i)& cdat< cvec(i+1)),'Facecolor',newmap(i,:))
end


contour(Topo.lon.data(iy,ix),Topo.lat.data(iy,ix),Topo.topo.data(iy,ix),[0 0],'LineColor',[0.6 0.5 0.4]);

