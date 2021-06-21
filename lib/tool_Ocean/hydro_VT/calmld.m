%
% function [mlde,mlpt,mlps,mlpd,...
%          dephmld_err,tpotmld_err,psalmld_err,sig0mld_err,...
%	  dephmld_qc,tpotmld_qc,psalmld_qc,sig0mld_qc]=calmld(sig01,deph1,tpot1,psal1,...
%                     errsig01,errpres1,errtpot1,errpsal1,...
%		     qcsig01,qcpres1,qctpot1,qcpsal1,critere)
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
%                        - sigm(psal(1er niv),tpot(1er niv)-0.02,0)
%
% le 1er niveau est egal a : 1er niv = 10 m
%
% creation V. Thierry juin 2005
%

function [dephmld,tpotmld,psalmld,sig0mld,dephmld_err,tpotmld_err,psalmld_err,sig0mld_err,dephmld_qc,tpotmld_qc,psalmld_qc,sig0mld_qc]=calmld(sig01,deph1,tpot1,psal1,errsig01,errpres1,errtpot1,errpsal1,qcsig01,qcpres1,qctpot1,qcpsal1,critere)

len=length(sig01);
ifin=find(finite(sig01)&deph1>=0);
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
  deph1(ifin(1))
end


if strcmp(casnan,'n') == 1
  dephmld=deph1(indmld);
  sig0mld=mynanmean(sig01(1:indmld));
  tpotmld=mynanmean(tpot1(1:indmld));
  psalmld=mynanmean(psal1(1:indmld));
  
  dephmld_err=errpres1(indmld);
  sig0mld_err=mynanmean(errsig01(1:indmld));
  tpotmld_err=mynanmean(errtpot1(1:indmld));
  psalmld_err=mynanmean(errpsal1(1:indmld));
  
  dephmld_qc=qcpres1(indmld);
  sig0mld_qc=max(qcsig01(1:indmld));
  tpotmld_qc=max(qctpot1(1:indmld));
  psalmld_qc=max(qcpsal1(1:indmld));
  
elseif strcmp(casnan,'y') == 1
  dephmld=NaN;dephmld_err=NaN;dephmld_qc=NaN;
  sig0mld=NaN;sig0mld_err=NaN;sig0mld_qc=NaN;
  tpotmld=NaN;tpotmld_err=NaN;tpotmld_qc=NaN;
  psalmld=NaN;psalmld_err=NaN;psalmld_qc=NaN;
end
