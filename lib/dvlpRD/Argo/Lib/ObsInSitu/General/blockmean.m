function [VARGD]=blockmean(lon,lat,VAL,domain)
% function [VARGD]=blockmean(lon,lat,VAL,domain)
% -----------------------------------
%   [OUTPUTS]=blockmean(INPUTS)
%   PURPOSE : compute the average of a xyz data on a regular grid, using the GMT  function blockmean
%   AUTHOR  : ccabanes
%   DATE    : 31/07/07
%   CALLED SUBROUTINES:
%
%   INPUT :
%%     VAL   
%             contains the scatter points
%      lon : longitude of scatter points
%      lat : latitude of scatter points
%
%%     domain (structure)
%             contains the domain info .latmin
%                                      .latmax
%                                      .lonmin
%                                      .lonmax
%                  and grid spacing    .spacex
%                                      .spacey
%
%   OUTPUT :
%%     VARA   (structure)
%             contains the averaged variable and grid description
% -------------------------------------
%
%
%% Simple average

ksel=find(~isnan(VAL)&~isnan(lon)&~isnan(lat));

%keyboard  
A(:,1)= lon(ksel);
A(:,2)= lat(ksel);
A(:,3)= VAL(ksel);
%A(:,4)= cos(A(:,2)*pi/180);
save ('/tmp/tempucor.dat','A','-ascii')

% block average
expre=['!blockmean /tmp/tempucor.dat -R' num2str(domain.lonmin) '/' num2str(domain.lonmax) '/' num2str(domain.latmin) '/' num2str(domain.latmax)];
expre=[expre ' -I' num2str(domain.spacex) '/' num2str(domain.spacey) ' -C  -Lx >/tmp/tempucor.int']
%expre=[expre ' -I' num2str(domain.spacex) '/' num2str(domain.spacey) '  -C -L -W >/tmp/tempucor.int']
eval(expre)
load ('/tmp/tempucor.int', '-ascii')
B=tempucor;

newlat=domain.latmin:domain.spacey:domain.latmax;
newlon=domain.lonmin:domain.spacex:domain.lonmax;

VARA.gd.dimy=length(newlat);
VARA.gd.dimx=length(newlon);

VARA.gd.lon=newlon;
VARA.gd.lat=newlat;

VARA.field=NaN*zeros(VARA.gd.dimy,VARA.gd.dimx);

for kk=1:length(B)
    jj=find(VARA.gd.lat==B(kk,2));
    ii=find(VARA.gd.lon==B(kk,1));
    VARA.field(jj,ii)=B(kk,3);
end

VARGD=VARA;