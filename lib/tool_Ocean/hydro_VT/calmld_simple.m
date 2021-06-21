%
% function
% [dephmld,tpotmld,psalmld,sig0mld]=calmld_simple(sig01,deph1,tpot1,psal1,critere)
%
% Function qui calcule la profondeur de la
% couche de melange en fonction de divers
% criteres
% critere = 
% 'dsig003'      : sig0(1er niv) - sig0(base mld) <= 0.03
% 'dt02'         : tpot(1er niv) - tpot(base mld) <= 0.2
% 'dsig003&dt02' : sig0(1er niv) - sig0(base mld) <= 0.03
%               et tpot(1er niv) - tpot(base mld) <= 0.2
% 'dsigdt02'     : sig0(1er niv) - sig0(base mld) <= dsigm
%               et tpot(1er niv) - tpot(base mld) <= 0.2
%             avec dsigm = sigm(psal(1er niv),tpot(1er niv),0) 
%                        - sigm(psal(1er niv),tpot(1er niv)-0.2,0)
%
% le 1er niveau est egal a : 1er niv = 10 m
%
% creation V. Thierry juin 2005
%

function [dephmld,tpotmld,psalmld,sig0mld]=calmld_simple(sig01,deph1,tpot1,psal1,critere)

len=length(sig01);
ifin=find((isfinite(sig01)==1) & (deph1>=0));
casnan='n';

if isempty(ifin) == 1
  casnan='y';
elseif deph1(ifin(1))<=10 % 1er niveau fini doit etre == 10 m
  iz=ifin(1);
  indmld=iz;
  sig20=sig01(iz);
  tpo20=tpot1(iz);
  psa20=psal1(iz);
  conti='y';
  while strcmp(conti,'y') == 1
    iz=iz+1;
    indmld=iz-1;
    if iz <= len
      
      switch critere
	
	case 'dsig003'
	  signp1=sig01(iz);
          if isfinite(signp1) & (signp1 <= sig20+0.03)
            conti='y';
          else
	    conti='n';
	  end
	  
	case 'dt02'
	  tponp1=tpot1(iz);
	  if isfinite(tponp1) & (abs(tponp1-tpo20)<= 0.2)
            conti='y';
          else
	    conti='n';
	  end
	  
	case 'dsig003&dt02'
	  signp1=sig01(iz);
	  tponp1=tpot1(iz);
	  
	  if isfinite(tponp1) & (abs(tponp1-tpo20) <= 0.2) & (signp1 <= sig20+0.03)
            conti='y';
          else
	    conti='n';
	  end
	  
	case 'dsigdt02'
	  [svan,sigma1]=swstat90(psa20,tpo20,0);
	  [svan,sigma2]=swstat90(psa20,tpo20-0.2,0);
	  dsigm=abs(sigma2-sigma1);
	  signp1=sig01(iz);
	  tponp1=tpot1(iz);
	  
	  if isfinite(tponp1) & (abs(tponp1-tpo20) <= 0.2) & (signp1 <= sig20+dsigm)
            conti='y';
          else
	    conti='n';
	  end
	  
      end
 
    else
      conti='n';
    end
  end
else
  casnan='y';
  deph1(ifin(1));
end


if strcmp(casnan,'n') == 1
  dephmld=deph1(indmld);
  sig0mld=mynanmean(sig01(1:indmld));
  tpotmld=mynanmean(tpot1(1:indmld));
  psalmld=mynanmean(psal1(1:indmld));
  
  
elseif strcmp(casnan,'y') == 1
  dephmld=NaN;
  sig0mld=NaN;
  tpotmld=NaN;
  psalmld=NaN;
end
