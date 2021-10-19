%  --------------------------------------------------------------------------
% function [lat_y,lon_x] = ...
%                  cylgeo(lat_orig, lon_orig, X_coor, Y_coor, l_axis, s_axis)
%
%  --------------------------------------------------------------------------
%  Conversion from cylindrical(x,y) to geographical (long,lat) coordinates
%
%   Input : 
%       lat_orig : origine of latitude for projection (deg.dec)
%       lon_orig : origine of longitude for projection (deg.dec)
%       X_coor   : coordinates relative to lat_orig and lon_orig (in meters) 
%       Y_coor   :                     "
%       l_axis   : half long axis for the reference ellipsoid
%       s_axis   : half short axis for the reference ellipsoid 
%
%   Output : 
%       lat_y 	: latitude (deg.dec)
%       lon_x	: longitude (deg.dec)

%
%  V 1.01     					A. Guillot,  	juin 1993
%  V 1.02     					F. Gaillard, 	Octobre 1998 
%  V 1.03                       		T. Terre, 	Mars 1999,
%          Change comments, add reference elliposid parameters
%  V 1.04 					F. Gaillard,	December 1999
%  	  Allows for vectors of coordinates	
%  --------------------------------------------------------------------------
%

function [lat_y,lon_x] = ...
                  cylgeo(lat_orig, lon_orig, X_coor, Y_coor, l_axis, s_axis)

%
% Set constant relative to ellipsoid and prepares correction coefficients

geocst

%  Computes P1
%  --------------
B  = (Y_coor - y0)/K0 + mo;
P0 = B/c1;
B0 = c1*P0 + c2*sin(2*P0) + c3*sin(4*P0);
epsil = 1.0;
while (epsil > 1.0e-12);
    P1 = P0 + (B-B0)/c1;
    B1 = c1*P1+c2*sin(2*P1)+c3*sin(4*P1);
    B0 = B1;
    P0 = P1;
    epsil = max(abs((B-B1)/c1));
end
cos_P1 = cos(P1);
sin_P1 = sin(P1);


J2 = K0*l_axis*ones(size(X_coor))./sqrt(1.0 - e2*sin_P1.*sin_P1);
J3 = J2.*J2;

T  = sin_P1./cos_P1;
T2 = T.*T;
Z  = cos_P1.*cos_P1*e2/(1-e2);
Z2 = Z.*Z;
T2Z = T2.*Z;
QQ1 = 45*(2 + T2 - T2Z) - 162*Z;
QQ2 = 61 + 107*Z + T2.*QQ1;
QQ3 = 2*Z-3*Z2 - 1;
QQ4 = 28*T2 + 24*T2.*T2 + 8*T2Z + 6*Z + 5;
QQ5 = 5 - 3*T2.*QQ3 - 3*(Z2 - 2*Z);
QQ6 = -1 - 2*T2 - Z;

V   = X_coor - x0;
V2  = V.*V;
VV1 = V2./(12*J3);
VV2 = QQ5 - V2.*QQ2./(30*J3);
VV4 = 2*J3.*(VV1.*VV2 - (1+Z));

J0 = T.*V2./VV4;
J4 = QQ6 + V2.*QQ4./(20*J3);
J5 = (1 + V2.*J4./(6*J3));
J1 = V.*J5./(J2.*cos_P1);

%  Coordonnees du point en degres decimaux:
lat_y = (P1 + J0)/degrad;
lon_x = lon_orig + J1/degrad;



