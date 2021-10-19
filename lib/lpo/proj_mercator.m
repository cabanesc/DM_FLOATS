function [proj] = proj_mercator(y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Transformation d'une latitude par projection de Mercator
% Ellipsoide WGS 84
%
% V1.0 Y.A. 02/01/96 : Creation
%
% function [proj] = proj_mercator(y)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WGS 84
% a = 6378137.00;
% b = 6356752.314;
% e =sqrt(a^2-b^2)/a;

e=0.0818191913108718138;

rad=pi/180;

esinteta=e*sin(y*rad);

proj=log(tan(pi/4.0+(y*rad)/2.0)*(((1-esinteta)/(1+esinteta))^(e/2)))/rad;

return
