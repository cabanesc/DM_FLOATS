% function [vrpsig]=vorpotsigma(deph,pres,temp,psal,presref,tabsig,delta)
% 
% objet : Calcul de la vorticite potentielle sur des niveaux sigma
% -----
%
% entree : 
% ------
%          deph(nsta,ndep) : tableau des profondeurs (en m, positives)
%          pres(nsta,ndep) : tableau des pressions
%          temp(nsta,ndep) : tableau des temperatures in situ
%          psal(nsta,ndep) : tableau des salinite
%          lat(nsta)       : tableau des latitudes
%          presref : pression de reference 
%          tabsig  : vecteur des densites ou on souhaite
%          calculer la vorticite potentielle
%          delta   : difference de densite sur laquelle on va calculer
%          la vorticite, si delta est omis, delta =0.02
%
%
% sortie :
% ------
%          vrpsig  : tableau des vorticite potentielles calcules au
%          densite contenues dans tabsig
%
%

function [vrpsig]=vorpotsigma(deph,pres,temp,psal,presref,tabsig,lat,delta)

if  nargin < 8
  delta = 0.02;
end

fcor  = pi*sin(pi*lat/180.0)/(6*3600);

tabsig=tabsig(:);
[nsta,ndep]=size(pres);
lts=length(tabsig);
vrpsig=NaN*ones(nsta,lts);

[tpot0] = tetai(pres, temp, psal, presref);
[svan,sig02]=swstat90(psal,tpot0,presref);
st=size(tpot0);
if st(1) ~= nsta
  tpot0=tpot0';
  sig02=sig02';
end

sig02=correct_invert(sig02);

for ista=1:nsta
  inan=find(finite(sig02(ista,:)) & finite(deph(ista,:)));
  if isempty(inan) == 0 & length(inan) > 1
    dephsigp=interp1(sig02(ista,inan),deph(ista,inan)',(tabsig+delta/2)');
    dephsigm=interp1(sig02(ista,inan),deph(ista,inan)',(tabsig-delta/2)');
    vrpsig(ista,:)=[(fcor(ista)./(1000+tabsig)').*(delta./(dephsigp-dephsigm))];
  end
end
