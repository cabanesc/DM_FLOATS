%
% SODANO     Compute ranges on reference ellipsoid 
%
% CALL       [DIST] = SODANO(RLON1, RLAT1, RLON2, RLAT2, L_AXIS, S_AXIS)
%
% INPUT      RLON1    longitudes of starting points (deg.dec)
%            RLAT1    latitudes of starting points (deg.dec)
%            RLON2    longitudes of arrival points (deg.dec)
%            RLAT2    latitudes of arrival points (deg.dec)
%            L_AXIS   long axis for the reference ellipsoid
%            S_AXIS   short axis for the reference ellipsoid 
%
% OUTPUT     DIST     Ranges between starting and ending points (meters)
%

% HISTORY    Creation      F. Gaillard   juillet 1993
%            Vectorisation T.Terre   16/08/95
%            Add ellipsoid parameters  T.Terre  8/03/99
% REFERENCE  Cree a partir d'un routine EPSHOM
%

function [dist] = sodano(rlon1d,rlat1d,rlon2d,rlat2d, l_axis, s_axis)

deuxpi  = 2.d0*pi;
convrad = pi/180;
rlon1 = rlon1d*convrad;
rlon2 = rlon2d*convrad;
rlat1 = rlat1d*convrad;
rlat2 = rlat2d*convrad;

% Transforme les entrees en vecteur colonne
if (size(rlon1,1) < size(rlon1,2))
	rlon1 = rlon1';
	rlat1 = rlat1';
	rlon2 = rlon2';
	rlat2 = rlat2';
end;



a = l_axis;
b = s_axis;

flat  = 1.0 - b/a;
flat2 = flat*flat;
f1=flat2*1.25;
f2=flat2*0.50;
f3=flat2*0.25;
f4=flat2*0.125;
f5=flat2*0.0625;
f6=flat2+flat;
f7=f6+1.0;
f8=f6*0.5;

dist = zeros(length(rlon1),1);

beta1  = atan((1.-flat)*sin(rlat1)./cos(rlat1));
sbeta1 = sin(beta1);
cbeta1 = cos(beta1);
beta2  = atan((1.-flat)*sin(rlat2)./cos(rlat2));
sbeta2 = sin(beta2);
cbeta2 = cos(beta2);
dlat = rlat1 - rlat2;
dlon=rlon1-rlon2;
adell=abs(dlon);

ii = find(adell >= pi);
adell(ii) = deuxpi-adell(ii); 
sidel=sin(adell);
codel=cos(adell);
a1=sbeta1.*sbeta2;
b1=cbeta1.*cbeta2;
cophi=a1+b1.*codel;
tmp = (sbeta2.*cbeta1-sbeta1.*cbeta2.*codel);
tmp0 = (sidel.*cbeta2);
siphi=sqrt(tmp0.*tmp0+tmp.*tmp);

ii = find(siphi > 0 | siphi < 0);
c=b1(ii).*sidel(ii)./siphi(ii);
em=1.-c.*c;
phi=asin(siphi);
jj = find(cophi < 0);
phi(jj) = pi - phi(jj);

phisq=phi(ii).*phi(ii);
csphi=1./siphi(ii);
ctphi=cophi(ii)./siphi(ii);
psyco=siphi(ii).*cophi(ii);

term1=f7.*phi(ii);
term2=a1(ii).*(f6*siphi(ii)-f2*phisq.*csphi);
term3=em.*(f2*phisq.*ctphi-f8*(phi(ii)+psyco));
term4=a1(ii).*a1(ii).*f2.*psyco;
term5=em.*em.*(f5*(phi(ii)+psyco)-f2*phisq.*ctphi-f4*psyco.*cophi(ii).*cophi(ii));
term6=a1(ii).*em*f2.*(phisq.*csphi+psyco.*cophi(ii));
dist(ii)=b.*(term1+term2+term3-term4+term5+term6);

