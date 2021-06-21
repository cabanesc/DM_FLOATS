% function ds = andoyer(long1,lat1, long2,lat2)
% Calcul de la distance en km entre deux points sur l'ellipsoide
% a l'aide de la formule d'Andoyer.
% longitudes et latitudes en degres
% long2, lat2 sont soit de meme dimension
% que long1 et lat1, soit des scalaires

function ds = andoyer(long1,lat1, long2,lat2)

s1=size(long1);
s2=size(long2);

%if prod(s2) ~= 1 & prod(s2) ~= prod(s1)
%  disp(['Il y a un probleme de dimension des parametres d entree']);
%end

a = 6378137.;    % grand axe 
b = 6356752.314; % petit axe 
f = (a-b)/a;     % applatissement 

d=zeros(s1);
ds=zeros(s1);


l1 = lat1*pi/180;
g1 = long1*pi/180;
l2 = lat2*pi/180;
g2 = long2*pi/180;


d = acos(sin(l1).*sin(l2) + cos(l1).*cos(l2).*cos(g2-g1));

[iz jz]=find(long1==long2 & lat1==lat2);
div=1-cos(d);
div(iz,jz)=9999;

ds = 0.001 * a * (d - f * 0.25 * ((d-3*sin(d)) .* ...
    ((sin(l1) + sin(l2)).*(sin(l1) + sin(l2)))./(1+cos(d)) + ...
    (d+3*sin(d)) .* ((sin(l1) - sin(l2)).*(sin(l1) - sin(l2)))./div));

ds(iz,jz)=0;
