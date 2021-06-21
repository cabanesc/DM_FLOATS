function contour_bathy(Topo,VEC_REG)
% -========================================================
%   USAGE : contour_bathy(Topo,VEC_REG)
%   PURPOSE : trace les contours 0, 1000 2000 3000 et 4000 de la bathy dans la region VEC_REG
% -----------------------------------
%   INPUT :
%     IN1   (class)  -comments-
%             additional description
%     IN2   (class)  -comments-
%
%   OPTIONNAL INPUT :
%    OPTION1  (class)  -comments-
% -----------------------------------
%   OUTPUT :
%     OUT1   (class)  -comments-
%             additional description
%     OUT2   (class)  -comments-
%             additional description
% -----------------------------------
%   HISTORY  : created (2009) ccabanes
%            : modified (yyyy) byxxx 
%   CALLED SUBROUTINES: none
% ========================================================

hold on
box on
grid on
xlabel('longitude');ylabel('latitude')
if isfield(Topo,'lon')
    Topo.lon.data=Topo.lon.data(1:end-1);
    Topo.topo.data=Topo.topo.data(:,1:end-1);
    Topo=shiftEW(Topo,'lon','grwch');
    ii=find(Topo.lat.data>=VEC_REG(3)&Topo.lat.data<=VEC_REG(4));
    jj=find(Topo.lon.data>=VEC_REG(1)&Topo.lon.data<=VEC_REG(2));
    Topo.lon.data=Topo.lon.data(jj);
    Topo.lat.data=Topo.lat.data(ii);
    Topo.topo.data=Topo.topo.data(ii,jj);
    [junk,isort]=sort(Topo.lon.data);
    contour(Topo.lon.data(isort),Topo.lat.data,Topo.topo.data(:,isort),[0 0],'LineColor',[0 0 0]); 
    contour(Topo.lon.data(isort),Topo.lat.data,Topo.topo.data(:,isort),[-1000 -1000],'LineColor',[0.5 0.5 0.5]); 
    contour(Topo.lon.data(isort),Topo.lat.data,Topo.topo.data(:,isort),[-2000 -2000],'LineColor',[0.8 0.8 0.8]); 
    contour(Topo.lon.data(isort),Topo.lat.data,Topo.topo.data(:,isort),[-3000 -3000],'LineColor',[0.8 0.8 0.8]); 
    contour(Topo.lon.data(isort),Topo.lat.data,Topo.topo.data(:,isort),[-4000 -4000],'LineColor',[0.8 0.8 0.8]); 
elseif isfield(Topo,'grd_lon')
    Topo=shiftEW(Topo,'grd_lon','grwch');
    contour(Topo.grd_lon.data,Topo.grd_lat.data,Topo.topo.data,[0 0],'LineColor',[0 0 0]); 
    contour(Topo.grd_lon.data,Topo.grd_lat.data,Topo.topo.data,[-1000 -1000],'LineColor',[0.5 0.5 0.5]); 
    contour(Topo.grd_lon.data,Topo.grd_lat.data,Topo.topo.data,[-2000 -2000],'LineColor',[0.8 0.8 0.8]); 
    contour(Topo.grd_lon.data,Topo.grd_lat.data,Topo.topo.data,[-3000 -3000],'LineColor',[0.8 0.8 0.8]); 
    contour(Topo.grd_lon.data,Topo.grd_lat.data,Topo.topo.data,[-4000 -4000],'LineColor',[0.8 0.8 0.8]); 
end
axis(VEC_REG)