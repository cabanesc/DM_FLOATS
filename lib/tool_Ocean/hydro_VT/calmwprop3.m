%
% function [pres_mini,pres_mini_qc,pres_mini_err,...
%           sig0_mini,sig0_mini_qc,sig0_mini_err,...
%           tpot_mini,tpot_mini_qc,tpot_mini_err,...
%           psal_mini,psal_mini_qc,psal_mini_err,...
%           pres_mean,pres_mean_qc,pres_mean_err,...
%           sig0_mean,sig0_mean_qc,sig0_mean_err,...
%           tpot_mean,tpot_mean_qc,tpot_mean_err,...
%           psal_mean,psal_mean_qc,psal_mean_err,...
%           epais,epais_qc,epais_err]=calmwprop(sig01,pres1,tpot1,psal1,vorp1,...
%                                               errsig01,errpres1,errtpot1,errpsal1,...
%                                               qcsig01,qcpres1,qctpot1,qcpsal1,critere)
% 
%
% Function qui calcule les proprietes des eaux modales
% en fonction de divers
% criteres
% critere = 
% 'limvrp6'      : vrp < 6 1e-11 m-1.s-1
% 'limvrp5'      : vrp < 5 1e-11 m-1.s-1
% 'limvrp4'      : vrp < 4 1e-11 m-1.s-1*
%

function [pres_mini,pres_mini_qc,pres_mini_err,sig0_mini,sig0_mini_qc,sig0_mini_err,tpot_mini,tpot_mini_qc,tpot_mini_err,psal_mini,psal_mini_qc,psal_mini_err,pres_mean,pres_mean_qc,pres_mean_err,sig0_mean,sig0_mean_qc,sig0_mean_err,tpot_mean,tpot_mean_qc,tpot_mean_err,psal_mean,psal_mean_qc,psal_mean_err,epais,epais_qc,epais_err]=calmwprop3(dephmld,sig01,pres1,tpot1,psal1,vorp1,errsig01,errpres1,errtpot1,errpsal1,qcsig01,qcpres1,qctpot1,qcpsal1,critere)


if dephmld == 0
  dephmld=100;
end


pres_mini=NaN*ones(1,3);
pres_mini_qc=NaN*ones(1,3);
pres_mini_err=NaN*ones(1,3);

sig0_mini=NaN*ones(1,3);
sig0_mini_qc=NaN*ones(1,3);
sig0_mini_err=NaN*ones(1,3);

tpot_mini=NaN*ones(1,3);
tpot_mini_qc=NaN*ones(1,3);
tpot_mini_err=NaN*ones(1,3);

psal_mini=NaN*ones(1,3);
psal_mini_qc=NaN*ones(1,3);
psal_mini_err=NaN*ones(1,3);

pres_mean=NaN*ones(1,3);
pres_mean_qc=NaN*ones(1,3);
pres_mean_err=NaN*ones(1,3);

sig0_mean=NaN*ones(1,3);
sig0_mean_qc=NaN*ones(1,3);
sig0_mean_err=NaN*ones(1,3);

tpot_mean=NaN*ones(1,3);
tpot_mean_qc=NaN*ones(1,3);
tpot_mean_err=NaN*ones(1,3);

psal_mean=NaN*ones(1,3);
psal_mean_qc=NaN*ones(1,3);
psal_mean_err=NaN*ones(1,3);

epais=NaN*ones(1,3);
epais_qc=NaN*ones(1,3);
epais_err=NaN*ones(1,3);

sig0_std=NaN*ones(1,3);
sig0_levmin=NaN*ones(1,3);
sig0_levmax=NaN*ones(1,3);
pres_levmin=NaN*ones(1,3);
pres_levmax=NaN*ones(1,3);

  
% Definition des limites

limvrp6=6*1e-11;
limvrp5=5*1e-11;
limvrp4=4*1e-11;
si0bound=[26.85 27.7];
prsbnd=[dephmld 900];
salbnd=[35 36];

% Critere
switch critere
  case 'limvrp6'
    isout_low_vrp=find(vorp1>limvrp6 | ...
	               sig01>=si0bound(2) | sig01<=si0bound(1) | ...
		       pres1<prsbnd(1) | ...
		       psal1<salbnd(1) | psal1>salbnd(2));
		   %| pres1>prsbnd(2)
  case 'limvrp5'
    isout_low_vrp=find(vorp1>limvrp5 | ...
	               sig01>=si0bound(2) | sig01<=si0bound(1) | ...
		       pres1<prsbnd(1) | ...
		       psal1<salbnd(1) | psal1>salbnd(2));
		   %| pres1>prsbnd(2)
  case 'limvrp4'
    isout_low_vrp=find(vorp1>limvrp4 | ...
	               sig01>=si0bound(2) | sig01<=si0bound(1) | ...
		       pres1<prsbnd(1) | ...
		       psal1<salbnd(1) | psal1>salbnd(2));
		   %| pres1>prsbnd(2)
end
    
vorpini=vorp1;
vorpini(isout_low_vrp)=NaN;
tablevini=sig01;
isnofin=find(isnan(vorpini));
tablevini(isnofin)=NaN;
isfin=find(isfinite(vorpini));

if isempty(isfin) == 0
  
  difin=find(diff(isfin)>1);
  st=size(difin);
  if st(1)>st(2)
    difin=difin';
  end
  
  tabideb=[1 difin+1];
  tabifin=[difin length(isfin)];
  thi=tabifin-tabideb+1;
  [vsort,isort]=sort(thi,'descend');
  tabiimax=isort(1:min(3,length(tabifin)));
  
  for iimax=1:min(3,length(tabifin))
    vorp2=NaN*ones(size(vorpini));
    sig02=NaN*ones(size(tablevini));
    ideb=tabideb(tabiimax(iimax));
    ifin=tabifin(tabiimax(iimax));
    vorp2(isfin(ideb:ifin))=vorpini(isfin(ideb:ifin));
    sig02(isfin(ideb:ifin))=tablevini(isfin(ideb:ifin));
  
    % minimum de VORP
    [valmini levvrpmini]=min(vorp2);
    % pression, densite, tpot, psal ou VORP est minimum
    pres_mini(iimax)=pres1(levvrpmini);
    sig0_mini(iimax)=sig01(levvrpmini);
    tpot_mini(iimax)=tpot1(levvrpmini);
    psal_mini(iimax)=psal1(levvrpmini);
  
    % densite min et max ou VORP est inferieur a une valeur limite
    sig0_levmin(iimax)=sig02(isfin(ideb));
    sig0_levmax(iimax)=sig02(isfin(ifin));

    % pression min et max ou VORP est inferieur a une valeur limite
    % Epaisseur
    deph_levmin(iimax)=pres1(isfin(ideb));
    deph_levmax(iimax)=pres1(isfin(ifin));
    epais(iimax)=deph_levmax(iimax)-deph_levmin(iimax);
  
    % densite moyenne et std ou VORP est inferieur a une valeur limite
    sig0_mean(iimax)=mean(sig02(isfin(ideb:ifin)));
    sig0_std(iimax)=std(sig02(isfin(ideb:ifin)));
  
    % tpot et psal moyenne et std ou VORP est inferieur a une valeur limite
    tpot_mean(iimax)=mean(tpot1(isfin(ideb:ifin)));
    psal_mean(iimax)=mean(psal1(isfin(ideb:ifin)));
    pres_mean(iimax)=mean(pres1(isfin(ideb:ifin)));

  
    % QC et erreur
    pres_mini_qc(iimax)=qcpres1(levvrpmini);
    pres_mini_err(iimax)=errpres1(levvrpmini);
  
    sig0_mini_qc(iimax)=qcsig01(levvrpmini);
    sig0_mini_err(iimax)=errsig01(levvrpmini);
  
    tpot_mini_qc(iimax)=qctpot1(levvrpmini);
    tpot_mini_err(iimax)=errtpot1(levvrpmini);
  
    psal_mini_qc(iimax)=qcpsal1(levvrpmini);
    psal_mini_err(iimax)=errpsal1(levvrpmini);
  
    epais_qc(iimax)=qcpres1(levvrpmini);
    epais_err(iimax)=errpres1(levvrpmini);
    
    pres_mean_qc(iimax)=qcpres1(levvrpmini);
    pres_mean_err(iimax)=errpres1(levvrpmini);
  
    sig0_mean_qc(iimax)=qcsig01(levvrpmini);
    sig0_mean_err(iimax)=errsig01(levvrpmini);
  
    tpot_mean_qc(iimax)=qctpot1(levvrpmini);
    tpot_mean_err(iimax)=errtpot1(levvrpmini);
  
    psal_mean_qc(iimax)=qcpsal1(levvrpmini);
    psal_mean_err(iimax)=errpsal1(levvrpmini);
  end 
end

  
  