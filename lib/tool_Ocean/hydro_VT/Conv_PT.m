%-------------------------------------------------------------------------------
%
% Conv_PT        - Calcule pressure and in-situ temperature from
%		   Z, Tpot, S 
%
%  function[Pres, Temp, Sigi] = Conv_PT(Zzz, Tpot, Sali, Ylat);
%
%-------------------------------------------------------------------------------
%                                                  
%     input : 
%     ------   
%        Zzz 	: vector of depth (in meters) 
%        Tpot 	: vector of potential temperature (in degrees C)
%        Sal 	: vector of salinity (in PSU)
%        Ylat 	: Latitude in decimal degrees
%  
%     output : 
%     ------   
%        Pres     : vector of presure (in db)
%        Temp     : vector of in-situ temperature (in degrees C)
%        Sigi     : vector of in-situ density anomaly (in kg m**-3)
%
%     internal calls to subroutines : 
%     -----------------------------  
%     tetai, swstat90, zenprs 
%
% Version:
% -------
%  1.01 Création (d'après clip_bsv_moy)          14/03/2001 F. Gaillard
% 
%-------------------------------------------------------------------------------

function[Pres, Temp, Sigi] = Conv_PT(Zzz, Tpot, Sali, Ylat);

pref = zeros(size(Tpot));			
pres0 = abs(Zzz);			% Premieres approximations
zcoo = -pres0;					  
temp0 = tetai(pref, Tpot, Sali, pres0);	
[xbid, sigma0] = swstat90 (Sali, temp0, pres0);

epsilon = 1;

while epsilon > 1e-4			%  Iteration
   [pres1] = zenprs(zcoo, sigma0, Ylat);
   temp1 = tetai(pref, Tpot, Sali, pres1);	
   [xbid, sigma1] = swstat90 (Sali, temp1, pres1);
   epsilon = max(abs(sigma1-sigma0));
   sigma0 = sigma1;
end

Sigi = sigma1;
Pres = pres1;
Temp = temp1;
