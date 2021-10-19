%  GEOCST           Compute projection constants based on the reference
%                    ellipsoid. L_AXIS and S_AXIS are respectively the 
%                    long axis and short axis of the reference ellipsoid. 
%
%

% V1.01             F.Gaillard
% V1.02  8/03/1999  T.Terre  Set the reference ellipsoid characteristics as
%                            parameters

 
degrad = pi/180;

%  ellipsoide de référence = WGS 84
%  --------------------------------
%l_axis = 6378137.00;	
%s_axis = 6356752.314;

flat  = 1.0 - s_axis/l_axis;
FM1 = 1.0/flat;		% inverse aplatissement 1/f

%  projection cylindrique : 
%  ----------------------
K0 = 1.0;	% coefficient réducteur d'échelle au méridien origine
x0 = 0;		% coordonnées x,y à l'origine de projection
y0 = 0;
L0 = lat_orig*degrad;


%*************************************************************************/
% calcul des constantes de la projection cylindrique pour une ellipsoide */
% donnée et une projection donnée                                        */
%									  */
%			    e2  : carré de l'excentricité		  */
%			    c1  :					  */
%			    c2  : > constantes pour le développement      */
%			    c3  :   du méridien				  */
%			    mo  : dévlpt méridien origine projection 	  */
%*************************************************************************/

e2 = 1/FM1*(2-1/FM1);
j0 = l_axis*(1- e2);
j1 = 3*e2;
c1 = j0*(1+j1/4*(1+(15*e2)/16*(1+(35*e2)/36)));
c2 =-j0*(1+(5*e2)/4*(1+(35*e2)/32))*j1/8;
c3 = j0*(1+(7*e2)/4)*(15*e2*e2)/256;
mo = c1*L0+c2*sin(2*L0)+c3*sin(4*L0);
