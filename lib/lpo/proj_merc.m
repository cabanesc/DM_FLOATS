function [proj] = proj_mercator(y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Transformation d'une latitude par projection de Mercator
% Ellipsoide WGS 84
%
% V1.0 Y.A. 02/01/96 : Creation
% V1.1 F.G. 09/10/96 : Modification, traitement des vecteurs
%
% function [proj] = proj_mercator(y)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WGS 84
% a = 6378137.00;
% b = 6356752.314;
% e =sqrt(a^2-b^2)/a;

e = 0.0818191913108718138;

rad  = pi/180;

yrad     = y*rad;
esinteta = e*sin(yrad);
rapsin   = ((1-esinteta)./(1+esinteta));
rapsexp  = rapsin.^(e/2);
tanyrad  = tan(pi/4.0 + yrad/2.0);

proj     = log(tanyrad.*rapsexp)/rad;
