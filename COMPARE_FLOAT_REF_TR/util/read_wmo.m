function [CTD] = read_wmo(wmo_file,ctd_to_exclude,base);

if strcmp(lower(base),'argo')
p=strsplit(wmo_file,'argo_') ;
num_wmo=str2num(p{3});
elseif strcmp(lower(base),'ctd')
p=strsplit(wmo_file,'ctd_') ;
num_wmo=str2num(p{2});
end
wmo=load(wmo_file);
  %keyboard  
to_use=[];
%to_use = find(ismember(strtok(strtrim(wmo.source),'_'),strtrim(ctd_to_exclude))==0);
%not_use = find(ismember(strtok(strtrim(wmo.source),'_'),strtrim(ctd_to_exclude))==1);
to_use = find(ismember(strtrim(wmo.source),strtrim(ctd_to_exclude))==0);
not_use = find(ismember(strtrim(wmo.source),strtrim(ctd_to_exclude))==1);
%not_use2 = find(ismember(strtrim(wmo.source),strtrim(ctd_to_exclude))==1);

wmo.lat=wmo.lat(to_use);
wmo.long=wmo.long(to_use);
wmo.dates=wmo.dates(to_use);
wmo.source=wmo.source(to_use);
wmo.sal=wmo.sal(:,to_use);
wmo.pres=wmo.pres(:,to_use);
wmo.temp=wmo.temp(:,to_use);
wmo.ptmp=wmo.ptmp(:,to_use);

ctddatesstr=num2str(wmo.dates');
wmo.dates=datenum(str2num(ctddatesstr(:,1:4)),str2num(ctddatesstr(:,5:6)),str2num(ctddatesstr(:,7:8)));
wmo.profil=[1:length(wmo.lat)];
[date_sort,isort]=sort(wmo.dates);
%keyboard
nprf=length(wmo.lat);
CTD.pres.data = wmo.pres(:,isort)';
CTD.psal.data = wmo.sal(:,isort)';
CTD.temp.data = wmo.temp(:,isort)';
CTD.tpot.data = wmo.ptmp(:,isort)';
CTD.longitude.data = wmo.long(isort)';
CTD.latitude.data= wmo.lat(isort)';
CTD.dates.data = wmo.dates(isort);
CTD.boxes.data = num_wmo*ones(nprf,1);
CTD.source.data = wmo.source(isort);
%ctdsources = wmo.source;
[null,CTD.sig0.data]=swstat90(CTD.psal.data,CTD.tpot.data,0);
if size(CTD.pres.data) ~= size(CTD.sig0.data)
CTD.sig0.data= CTD.sig0.data';
end
%[CTD.temp2.data]= sw_temp(CTD.psal.data,CTD.tpot.data,CTD.pres.data,0);

CTD=shiftEW(CTD,'longitude','grwch');


% mise en forme des tableaux
%CTD.latitude.data = repmat(CTD.latitude.data,1,size(CTD.psal.data,2));
%CTD.longitude.data = repmat(CTD.longitude.data,1,size(CTD.psal.data,2)); 
CTD.profil.data = repmat([1:size(CTD.psal.data,1)]',1,size(CTD.psal.data,2));
CTD.profil_orig.data = wmo.profil(isort)';
%CTD.dates.data = repmat(CTD.dates.data,1,size(CTD.psal.data,2));
%CTD.boxes.data = repmat(CTD.boxes.data,1,size(CTD.psal.data,2));