% function [yearind,naoind]=getnao(type)
% Fonction qui lit l'indice NAO
% entree
%       type = 'hurrelldta', 'hurrelpc' ou 'ncep'
%       si type = '' alors type = 'hurrellpc'
%
% sortie
%       yearind = annee
%       naoind = index NAO
%
% l'indice NAO de HURREL est calculé pour les mois DJFM
% l'annee correspondante est celle de JFM
%
% l'indice NAO de NCEP est la moyenne des mois de DJFM
% 
% voir /home/revellata/vthierry/ATLNORD/NAO/
%

function [yearind,naoind]=getnao(type)

inpath='/home1/armen/perso/vthierry/matlab/NAO/';

if isempty(type) == 1
  type='hurrel';
end

if strcmp(type,'hurreldta') == 1
  load([inpath 'naodjfmindex.1864-2005.asc']);
  yearind=naodjfmindex(:,1);
  naoind=naodjfmindex(:,2);
elseif strcmp(type,'hurrelpc') == 1
  tabnao=load([inpath 'naodjfm_pcbased_1899-2008_hurrel.asc']);  
  [nl,nc]=size(tabnao);
  ind=0;
  for iline=2:nl
    for icol=2:nc
      ind=ind+1;
      yearind(ind)=tabnao(iline,1)+tabnao(1,icol);
      naoind(ind)=tabnao(iline,icol);
    end
  end
  naoind(find(naoind==-99))=NaN;
  inonan=find(isfinite(naoind));
  naoind=naoind(inonan);
  yearind=yearind(inonan);
elseif strcmp(type,'ncep') == 1
  tabdat=load([inpath 'ncep_norm.nao.monthly.b5001.current.ascii']);
  tabyear=tabdat(:,1);
  tabmonth=tabdat(:,2);
  tabnao=tabdat(:,3);
  indy=0;
  for iyear=tabyear(1):tabyear(end)
    indy=indy+1;
    if iyear == tabyear(1)
      indmean=find(tabyear==iyear & tabmonth<=3);
    else
      indmean=find((tabyear==iyear & tabmonth<=3) | (tabyear==iyear-1 & tabmonth==12));
    end
    naoind(indy)=mean(tabnao(indmean));
    yearind(indy)=iyear;
  end
else
  disp('type non valable');
end
