% [delta_upbot,delta_botup,is_inv_dens,rhop_i,rhop_ip1,Pref] = density_test(P,T,S)
%
% P(1,nlevel) pression
% T(1,nlevel) T in situ 
% S(1,nlevel) S in situ
%
% delta_upbot (1,nlevel) : test up to bottom : si delta_upbot >= 0.03kg/m3 flag 4 du niveau i
% delta_botup (1,nlevel) : test bottom to up : delta_botup <= -0.03kg/m3 flag 4 du niveau i
% is_inv_dens (1,nlevel) : is_inv_dens=1 si flag 4 du niveau i; ==0 sinon

function [delta_upbot,delta_botup,is_inv_dens,rhop_i,rhop_ip1,Pref] =test_densite14(P,T,S)

% Addpath commented by T. Reynaud 07/09/2020
%addpath('/home1/triagoz/matlab/outils_matlab/seawater/seawater_330_its90_lpo');

if size(P,1)> 1
  P=P';
end
if size(T,1)> 1
  T=T';
end
if size(S,1)> 1
  S=S';
end


% De haut en bas controle entre deux niveaux consecutifs i et i+1 (Pi < Pi+1 )

Pref = (P(1:end-1)+P(2:end))/2;


%calul de rhoi (Pi, Ti, Si, Pref) referencee au niveau Pref 

rhop_i = sw_pden(S(1:end-1),T(1:end-1),P(1:end-1),Pref);

% calul de rhoi+1 (Pi+1, Ti+1, Si+1, Pref) referencee au niveau Pref 

rhop_ip1 = sw_pden(S(2:end),T(2:end),P(2:end),Pref);

delta_upbot = [rhop_i-rhop_ip1, NaN]; % si delta_upbot >= 0.03kg/m3 flag 4 du niveau i


% De bas en haut controle entre deux niveaux consecutifs i et i+1 (Pi > Pi+1 )


Pref = (P(2:end)+P(1:end-1))/2;


%calcul de rhoi (Pi, Ti, Si, Pref) referencee au niveau Pref 

rhop_i = sw_pden(S(2:end),T(2:end),P(2:end),Pref);

% calcul de rhoi+1 (Pi+1, Ti+1, Si+1, Pref) referencee au niveau Pref 

rhop_ip1 = sw_pden(S(1:end-1),T(1:end-1),P(1:end-1),Pref);

delta_botup = [NaN, rhop_i-rhop_ip1]; % si delta_botup <= -0.03kg/m3 flag du niveau i

is_inv_dens = delta_upbot>=0.03 | delta_botup<=-0.03;


