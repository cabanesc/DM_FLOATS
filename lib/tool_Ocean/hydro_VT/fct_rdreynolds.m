% function [tablon,tablat,tabdate,tabsst,tabice]=fct_rdreynolds(ficname)
% Fonction qu lit la SST Reynolds
% les donnees sont sous
% /home1/revellata/vthierry/REYNOLDS/
% 

function [tablon,tablat,tabdate,tabsst,tabice]=fct_rdreynolds(ficname)

inpath='/home1/revellata/vthierry/REYNOLDS/';

fid=fopen([inpath 'lstags.onedeg.dat'],'r','b');
meshmask=reshape(fread(fid,360*180,'float32'),360,180);
fclose(fid);

tablon=[0.5:359.5];
tablat=-89.5:89.5;

fid=fopen([inpath ficname],'r','b');

if strcmp('oiv2mon.2008',ficname) == 1
  nmonth=2;
else
  nmonth=12;
end

for imonth=1:nmonth
  % record 1
  nrec=fread(fid,1,'int32');
  tabdate(:,imonth)=fread(fid,7,'int32');
  index(imonth)=fread(fid,1,'int32');
  nrec=fread(fid,1,'int32');
  % record 2
  nrec=fread(fid,1,'int32');
  clear tmp1
  tmp1=reshape(fread(fid,360*180,'float32'),360,180);
  tmp1(find(meshmask==0))=NaN;
  tabsst(:,:,imonth)=tmp1;
  nrec=fread(fid,1,'int32');
  % record 3
  nrec=fread(fid,1,'int32');
  clear tmp2
  tmp2=reshape(fread(fid,360*180,'char'),360,180);
  tmp2(find(meshmask==0))=NaN;
  tabice(:,:,imonth)=tmp2;
  nrec=fread(fid,1,'int32');
end
fclose(fid);
