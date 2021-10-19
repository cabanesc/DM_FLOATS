%
% [X_COOR, Y_COOR] = GEOCYL(LAT_ORIG, LON_ORIG, LAT_Y, LON_X, L_AXIS, S_AXIS)
%  Conversion from geographical (long,lat) to cylindrical(x,y) coordinates
%
%   Input : 
%       lat_orig : origine of latitude for projection (deg.dec)
%       lon_orig : origine of longitude for projection (deg.dec)
%       lat_y 	 : latitude (deg.dec)
%       lon_x	 : longitude (deg.dec)
%       l_axis   : long axis for the reference ellipsoid
%       s_axis   : short axis for the reference ellipsoid 
%
%   Output : 
%       x_coor   : coordinates relative to lat_orig and lon_orig (in meters) 
%       y_coor   :                     "

%  V 1.01  Creation -- A. Guillot,  juin 1993
%  V 1.02  Update -- F. Gaillard, Octobre 1998
%  V 1.03  Change comments, add reference ellipsoid parameters -- T.Terre 03/99
%                  

function [x_coor, y_coor] = geocyl(lat_orig, lon_orig, lat_y, lon_x, l_axis, s_axis)

%
% Set constant relative to ellipsoid
%
geocst


lat_yr = lat_y*degrad;
b = c1*lat_yr + c2*sin(2*lat_yr) + c3*sin(4*lat_yr);
w = (lon_x - lon_orig)*degrad;
cos_l  = cos(lat_yr);
sin_l  = sin(lat_yr);
t      = sin_l/cos_l;
cos_l2 = cos_l*cos_l;
z      = e2/(1-e2)*cos_l2;

n=l_axis/sqrt((1-e2*sin_l*sin_l));
j0=w*w*cos_l2;
t2=t*t;
h = (n*j0*t)/2*(1+j0/12*(5-t2+z*(9+4*z)+j0/30+(61+t2*(t2-58)-z*(330*t2-270))));
p = 1+j0/6*(1-t2+z+j0/20*(5+t2*(t2-18)+z*(14-58*t2)));

y_coor = y0 + K0*(b-mo+h);
x_coor = x0 + K0*n*w*cos_l*p;

