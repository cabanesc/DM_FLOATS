% -========================================================
%   USAGE :   [msg,hf,h1,h2,h3]=cmp_prf_argo_ref(floatname,dacname,profnum,PARAM,CONFIG)
%   PURPOSE : plot data from a given float
% -----------------------------------
%   INPUT :
%    floatname  (char)  e.g. '690258'
%    dacname    (char)  e.g.  'coriolis'
%    numcycle   (float array)  e.g. 1 (first cycle 001) see option DIRECTION to consider descending profiles
%
%    PARAM   (structure)   input from verif_flag_1profil.m
%                           PARAM.REXTEND
%
%   CONFIG   (structure)  input from verif_flag_1profil.m
% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  :  created: V. Thierry - N. David - Juillet 2007
%               revised: C. Lagadec
%               revised: ccabanes - 2016
%   CALLED SUBROUTINES:
% ========================================================

function [msg,hf,h1,h2,h3]=cmp_prf_argo_ref(wmonum,dacname,profnum,PARAM,CONFIG)

if PARAM.CHECK_REF==1
    fw=fopen([CONFIG.DIR_PLOT '/log/log_bad_ref_' PARAM.REFERENCE '_' CONFIG.VERSION_OWC '.txt'],'a');
end


%if strcmp(PARAM.DATATYPE,'adj') == 1
FLOAT_SOURCE_NETCDF = CONFIG.(PARAM.DATAREP);
%elseif strcmp(PARAM.DATATYPE,'raw') == 1
%FLOAT_SOURCE_NETCDF = CONFIG.DIR_FTP;
%end
floatname = num2str(wmonum);

maxptvert=4100;
if strcmp(PARAM.REFERENCE,'ctd') == 1
    BOX=2;
elseif strcmp(PARAM.REFERENCE,'argo') == 1
    BOX=4;
end
NPROFMAX=100;

% Initialisation
argopath = FLOAT_SOURCE_NETCDF;

pathwmobox = [CONFIG.DIR_OWC CONFIG.VERSION_OWC '/data/'];

inpath=[argopath dacname '/' floatname '/profiles/'];

%float we want to exclude
ctd_to_exclude={};


% read the float data

pfnu=sprintf('%3.3i',profnum);
if strcmp(PARAM.DIRECTION,'D')
    pfnu=[pfnu 'D'];
end

fname=['D' floatname '_' pfnu '.nc'];
if exist([inpath fname]) == 2
    file_name=[inpath fname];
else
    fname=['R' floatname '_' pfnu '.nc'];
    if exist([inpath fname]) == 2
        file_name=[inpath fname];
    else
        disp(['Files ' fname  ' does not exist']);
        disp([inpath]);
        msg=['Profile does not exist'];
        h1=NaN;h2=NaN;h3=NaN;msg=NaN;hf=NaN;
        return
    end
end

disp(file_name)
[FL,DimL] = read_netcdf_allthefile(file_name);
FL= replace_fill_bynan(FL);
FL= format_flags_char2num(FL);

format_version=str2num(FL.format_version.data');
if format_version>=3
    is_primary=findstr_tab(FL.vertical_sampling_scheme.data,'Primary sampling');
    n_prof = find(is_primary);
    if isempty(n_prof)
        warning([file_in ' :Primary profile not found, n_prof=1 is considered'])
        n_prof=1;
    end
else
    n_prof=1;
end
[F,DIM] = extract_profile_dim(FL,DimL,'N_PROF',n_prof);

% Computation of useful varaibles
if strcmp(PARAM.DATATYPE,'adj')==1 % replace raw varaibles by adjusted variables
    F.pres=F.pres_adjusted;
    F.pres_qc=F.pres_adjusted_qc;
    F.temp_qc=F.temp_adjusted_qc;
    F.temp=F.temp_adjusted;
    F.psal=F.psal_adjusted;
    F.psal_qc=F.psal_adjusted_qc;
    
end


[n,m]=size(F.pres.data);
pr=reshape(F.pres.data,n*m,1);
te=reshape(F.temp.data,n*m,1);
ps=reshape(F.psal.data,n*m,1);
tp=tetai(pr(:),te(:),ps(:),0);
F.tpot.data=reshape(tp,n,m);
[null,tabsig0]=swstat90(ps,tp,0);
F.sig0.data=reshape(tabsig0,n,m);

F.tpot_qc.data=F.temp_qc.data;
F.sig0_qc.data=max([F.temp_qc.data;F.psal_qc.data]);

[F,DIM]=fill_structure(F.tpot.data,'TPOT',{'N_PROF','N_LEVELS'},F,DIM);
[F,DIM]=fill_structure(F.tpot_qc.data,'TPOT_QC',{'N_PROF','N_LEVELS'},F,DIM);
[F,DIM]=fill_structure(F.sig0.data,'SIG0',{'N_PROF','N_LEVELS'},F,DIM);
[F,DIM]=fill_structure(F.sig0_qc.data,'SIG0_QC',{'N_PROF','N_LEVELS'},F,DIM);

clear tabsig0 pr te ps tp

% look for wmo boxes that contains reference data close to the Argo profile

if strcmp(PARAM.REFERENCE,'ctd') == 1
    load([pathwmobox 'constants/wmo_boxes_ctd.mat']);
elseif strcmp(PARAM.REFERENCE,'argo') == 1
    load([pathwmobox 'constants/wmo_boxes_argo.mat']);
end

F=shiftEW(F,'longitude','pacif');

la_wmo_num=find_25boxes(F.longitude.data,F.latitude.data,la_wmo_boxes);
igood=find(la_wmo_num(:,BOX) == 1);
la_wmo_num=la_wmo_num(igood,:);
[nbox,null]=size(la_wmo_num);

F=shiftEW(F,'longitude','grwch');

initvar=NaN*ones(maxptvert,nbox*NPROFMAX);
[CTD,DIM_CTD]=fill_structure(initvar,'PRES',{'N_LEVELS','N_PROF'});
[CTD,DIM_CTD]=fill_structure(initvar,'PSAL',{'N_LEVELS','N_PROF'},CTD,DIM_CTD);
[CTD,DIM_CTD]=fill_structure(initvar,'TPOT',{'N_LEVELS','N_PROF'},CTD,DIM_CTD);
[CTD,DIM_CTD]=fill_structure(initvar,'TEMP',{'N_LEVELS','N_PROF'},CTD,DIM_CTD);
initvar=NaN*ones(nbox*NPROFMAX,1);
[CTD,DIM_CTD]=fill_structure(initvar,'LONGITUDE',{'N_PROF'},CTD,DIM_CTD);
[CTD,DIM_CTD]=fill_structure(initvar,'LATITUDE',{'N_PROF'},CTD,DIM_CTD);
[CTD,DIM_CTD]=fill_structure(initvar,'DATES',{'N_PROF'},CTD,DIM_CTD);
[CTD,DIM_CTD]=fill_structure(initvar,'BOXES',{'N_PROF'},CTD,DIM_CTD);
[CTD,DIM_CTD]=fill_structure(initvar,'NUMPB',{'N_PROF'},CTD,DIM_CTD);

ctdsources={};

for ibox=1:nbox
    if strcmp(PARAM.REFERENCE,'ctd') == 1
        % Modified by T. Reynaud 10/09/2020
        filewmo=[pathwmobox 'climatology/historical_ctd/ctd_' num2str(la_wmo_num(ibox,1)) '.mat'];
        if ~exist(filewmo)
            filewmo=[CONFIG.DIR_OWC CONFIG.DIR_OWC_CTD 'ctd_' num2str(la_wmo_num(ibox,1)) '.mat'];
        end
        wmo=load(filewmo);
    elseif strcmp(PARAM.REFERENCE,'argo') == 1
        % Modified by T. Reynaud 10/09/2020
        filewmo=[pathwmobox 'climatology/argo_profiles/argo_' num2str(la_wmo_num(ibox,1)) '.mat'];
        if ~exist(filewmo)
            filewmo=[CONFIG.DIR_OWC CONFIG.DIR_OWC_ARGO 'argo_' num2str(la_wmo_num(ibox,1)) '.mat'];
        end
        wmo=load(filewmo);
    end
    to_use =[1:length(wmo.source)];
    not_use =[];
    wmo.lat=wmo.lat(to_use);
    wmo.long=wmo.long(to_use);
    wmo.dates=wmo.dates(to_use);
    wmo.source=wmo.source(to_use);
    wmo.sal=wmo.sal(:,to_use);
    wmo.pres=wmo.pres(:,to_use);
    wmo.temp=wmo.temp(:,to_use);
    wmo.ptmp=wmo.ptmp(:,to_use);
    wmo.numpb=to_use;
    
    
    %%% use correlation scale to select data in the box (NPROFMAX in each box)
    longitude_large=1.5;
    latitude_large=1;
    phi_large=0.1;
    PV_float=compute_PV(F.longitude.data,F.latitude.data);
    PV_hist=compute_PV(wmo.long,wmo.lat);
    
    la_grid_long1=wmo.long;
    gg=find(wmo.long>180);
    la_grid_long1(gg)=wmo.long(gg)-360;
    
    correlation_large = (la_grid_long1-F.longitude.data).^2./longitude_large.^2 + (wmo.lat-F.latitude.data).^2./latitude_large.^2 +...
        ( (PV_float-PV_hist)./sqrt( PV_float.^2+PV_hist.^2 )./phi_large ).^2 ;
    
    distance = andoyer(wmo.long,wmo.lat,F.longitude.data,F.latitude.data);
    if PARAM.USE_PV==0
        [d2,isort] = sort(distance);
    else
        [d2,isort] = sort(correlation_large);
    end
    nprf = min(NPROFMAX,length(isort));
    [np,null] = size(wmo.pres);
    
    if np>maxptvert
        np=maxptvert;
    end
    %keyboard
    ideb = NPROFMAX*(ibox-1)+1;
    ifin = NPROFMAX*(ibox-1)+nprf;
    
    CTD.pres.data(1:np,ideb:ifin) = wmo.pres(1:np,isort(1:nprf));
    CTD.psal.data(1:np,ideb:ifin) = wmo.sal(1:np,isort(1:nprf));
    CTD.tpot.data(1:np,ideb:ifin) = wmo.ptmp(1:np,isort(1:nprf));
	CTD.temp.data(1:np,ideb:ifin) = wmo.temp(1:np,isort(1:nprf));
    CTD.longitude.data(ideb:ifin) = wmo.long(isort(1:nprf));
    CTD.latitude.data(ideb:ifin) = wmo.lat(isort(1:nprf));
    CTD.dates.data(ideb:ifin) = wmo.dates(isort(1:nprf));
    CTD.boxes.data(ideb:ifin) = la_wmo_num(ibox,1)*ones(nprf,1);
    CTD.numpb.data(ideb:ifin)=wmo.numpb(isort(1:nprf));
    ctdsources(ideb:ifin) = wmo.source(isort(1:nprf));
end

% select first 50 (PARAM.REXTEND) neighboring profiles
[null,CTD.sig0.data]=swstat90(CTD.psal.data,CTD.tpot.data,0);
if size(CTD.pres.data) ~= size(CTD.sig0.data)
    CTD.sig0.data= CTD.sig0.data';
end
[CTD,DIM_CTD]=fill_structure(CTD.sig0.data,'SIG0',{'N_LEVELS','N_PROF'},CTD,DIM_CTD);

CTD=shiftEW(CTD,'longitude','grwch');

clear d2 isort
longitude_large=1.5;
latitude_large=1;
phi_large=0.1;
PV_float=compute_PV(F.longitude.data,F.latitude.data);
PV_hist=compute_PV(CTD.longitude.data,CTD.latitude.data);

la_grid_long1=CTD.longitude.data;
gg=find(CTD.longitude.data>180);
la_grid_long1(gg)=CTD.longitude.data(gg)-360;

correlation_large = (la_grid_long1-F.longitude.data).^2./longitude_large.^2 + (CTD.latitude.data-F.latitude.data).^2./latitude_large.^2 +...
    ( (PV_float-PV_hist)./sqrt( PV_float.^2+PV_hist.^2 )./phi_large ).^2 ;

distance=andoyer(F.longitude.data,F.latitude.data,CTD.longitude.data,CTD.latitude.data);
if PARAM.USE_PV==0
    [d2,isort]=sort(distance);
else
[d2,isort]=sort(correlation_large);
end

index_sort_dist=isort;
[CTD,DIM_CTD]=extract_profile_dim(CTD,DIM_CTD,'N_PROF',index_sort_dist(1:PARAM.REXTEND));


% controle qualite succint de la base de reference
ideC1=~isnan(CTD.psal.data);
stdbase1=std(CTD.psal.data(ideC1));
meanbase1=mean(CTD.psal.data(ideC1));
isout1=abs(CTD.psal.data-meanbase1)>10*stdbase1;

if PARAM.CHECK_REF==1
    % ouliers are NaN;
    % log des problemes dans un fichiers log
    findout=find(sum(isout1,2)>=1);
    if isempty (findout)==0
        disp('')
        disp('WARNING')
        disp([num2str(length(findout)) ' Salinity outlier(s) are found in the reference database'])
        disp(['see: ' CONFIG.DIR_PLOT '/log/log_bad_' PARAM.REFERENCE '_' CONFIG.VERSION_OWC '.txt'])
        disp('')
    end
    for ilog=1:length(findout)
        [puob,ipsal_bad_value]=max(abs(CTD.psal.data(findout(ilog),:)-meanbase1));
        fprintf(fw,'%s\n',[num2str(CTD.boxes.data(findout(ilog))) ', ' num2str(CTD.numpb.data(findout(ilog))) ', ' num2str(CTD.psal.data(findout(ilog),ipsal_bad_value))]);
    end
    CTD.psal.data(isout1)=NaN;
    
end
index_closest_dist=1;

%  ilon=find(ctdlon>=260);
%  ctdlon(ilon)=ctdlon(ilon)-360;

% dates in the same format
ctddatesstr=num2str(CTD.dates.data);
CTD.dates.data=datenum(str2num(ctddatesstr(:,1:4)),str2num(ctddatesstr(:,5:6)),str2num(ctddatesstr(:,7:8)));
F.dates.data=F.juld.data+datenum(1950,1,1,0,0,0);
F.dates.FillValue_=99999;
[vdate,index_closest_date]=min(abs(CTD.dates.data-F.dates.data));


maxlon=max(max(CTD.longitude.data),F.longitude.data)+2;
minlon=min(min(CTD.longitude.data),F.longitude.data)-2;
maxlat=max(max(CTD.latitude.data),F.latitude.data)+2;
minlat=min(min(CTD.latitude.data),F.latitude.data)-2;

zone_visu=[minlat maxlat minlon maxlon];
reso='LR';proj='mercator';

% Figure LOC
%-----------

hf=figure;

plot_hist_loc(zone_visu,reso,proj,CTD,F,index_closest_date);

CTD_LON=CTD.longitude.data;
CTD_LAT=CTD.latitude.data;


% FIGURE (P,T) (P,S) (P,SIG0) CLOSEST
%------------------------------------

h1=figure;
orient landscape


subplot(2,2,1)
hold on, grid on, box on


reshapef=10;
plot_prof_closest(CTD,F,'pres','temp',index_closest_dist,index_closest_date,PARAM,'reshapex',reshapef);

xlabel('In situ temperature')
ylabel('Pressure (db)');


subplot(2,2,3)
hold on, grid on, box on

reshapef=10;
plot_prof_closest(CTD,F,'pres','psal',index_closest_dist,index_closest_date,PARAM,'reshapex',reshapef);

xlabel('Salinity')
ylabel('Pressure (db)');


subplot(1,2,2)
hold on, grid on, box on

reshapef=10;
plot_prof_closest(CTD,F,'pres','sig0',index_closest_dist,index_closest_date,PARAM,'reshapex',reshapef);

xlabel('Potential density')
ylabel('Pressure (db)');

if strcmp(PARAM.DATATYPE,'raw') == 1
    hs=suptitle({[num2str(wmonum) ' - Cycle' pfnu ' - Raw - Date Argo profile ' datestr(F.dates.data,'dd-mmm-yyyy')];...
        ['Dates historicals profiles ' datestr(CTD.dates.data(index_closest_dist)) ' (magenta) and ' datestr(CTD.dates.data(index_closest_date)) ' (blue)']});
elseif  strcmp(PARAM.DATATYPE,'adj') == 1
    hs=suptitle({[num2str(wmonum) ' - Cycle' pfnu ' - Adj - Date Argo profile ' datestr(F.dates.data,'dd-mmm-yyyy')];...
        ['Dates historicals profiles ' datestr(CTD.dates.data(index_closest_dist)) ' (magenta) and ' datestr(CTD.dates.data(index_closest_date)) ' (blue)']});
end
set(hs,'fontsize',11);


% FIGURE (P,T) (P,S) (P,SIG0) ALL
%------------------------------------


h2=figure;
orient landscape

subplot(2,2,1)
hold on, grid on, box on
mindepth=0;
reshapef=10;
[xmin,xmax,ymin,ymax]=plot_prof_hist(CTD,F,'pres','temp',index_closest_dist,index_closest_date,PARAM,'mindepth',mindepth,'reshapex',reshapef);
xlabel('In situ temperature')
ylabel('Pressure (db)');
plot_prof_closest(CTD,F,'pres','temp',index_closest_dist,index_closest_date,PARAM,'mindepth',mindepth,'reshapex',reshapef);


%set(hs,'fontsize',11);
ymin;
ymax;
%set(gca,'ylim',[ymin ymax]);
%set(gca,'xlim',[xmin xmax]);
if ~isnan(ymin)&~isnan(ymax)
    set(gca,'ylim',[ymin ymax]);
end
if ~isnan(xmin)&~isnan(xmax)
    set(gca,'xlim',[xmin xmax]);
end

subplot(2,2,3)
hold on, grid on, box on
reshapef=10;
[xmin,xmax,ymin,ymax]=plot_prof_hist(CTD,F,'pres','psal',index_closest_dist,index_closest_date,PARAM,'mindepth',mindepth,'reshapex',reshapef);

xlabel('Salinity')
ylabel('Pressure (db)');
plot_prof_closest(CTD,F,'pres','psal',index_closest_dist,index_closest_date,PARAM,'mindepth',mindepth,'reshapex',reshapef);
ymin;
ymax;
if ~isnan(ymin)&~isnan(ymax)
    set(gca,'ylim',[ymin ymax]);
end
if ~isnan(xmin)&~isnan(xmax)
    set(gca,'xlim',[xmin xmax]);
end
%set(gca,'ylim',[ymin ymax]);
%set(gca,'xlim',[xmin xmax]);

mindepth=PARAM.DEPTH_ZOOM;
subplot(1,2,2)
hold on, grid on, box on
resahpef=10;
[xmin,xmax,ymin,ymax]=plot_prof_hist(CTD,F,'pres','sig0',index_closest_dist,index_closest_date,PARAM,'mindepth',mindepth,'reshapex',reshapef);
plot_prof_closest(CTD,F,'pres','sig0',index_closest_dist,index_closest_date,PARAM,'mindepth',mindepth,'reshapex',reshapef);
xlabel('Potential density')
ylabel('Pressure (db)');
ymin;
ymax;
if ~isnan(ymin)&~isnan(ymax)
    set(gca,'ylim',[ymin ymax]);
end
if ~isnan(xmin)&~isnan(xmax)
    set(gca,'xlim',[xmin xmax]);
end
% set(gca,'ylim',[ymin ymax]);
% set(gca,'xlim',[xmin xmax]);
theax=gca;

if strcmp(PARAM.DATATYPE,'raw') == 1
    hs=suptitle({[num2str(wmonum) ' - Cycle ' pfnu ' - Raw - Date Argo profile ' datestr(F.dates.data,'dd-mmm-yyyy')];...
        ['Dates historicals profiles ' datestr(CTD.dates.data(index_closest_dist)) ' (magenta) and ' datestr(CTD.dates.data(index_closest_date)) ' (blue)']});
elseif  strcmp(PARAM.DATATYPE,'adj') == 1
    hs=suptitle({[num2str(wmonum) ' - Cycle D-' pfnu ' - Adj - Date Argo profile ' datestr(F.dates.data,'dd-mmm-yyyy')];...
        ['Dates historicals profiles ' datestr(CTD.dates.data(index_closest_dist)) ' (magenta) and ' datestr(CTD.dates.data(index_closest_date)) ' (blue)']});
end
hs.FontWeight='bold';
colorbar(theax)

%FIGURE (T,S) CLOSEST and ALL
%------------------------------------

h3=figure;
orient landscape

subplot(4,8,[1:2,9:10])
hold on, grid on, box on
mindepth=0;
[xmin,xmax,ymin,ymax]=plot_prof_hist(CTD,F,'tpot','psal',index_closest_dist,index_closest_date,PARAM,'mindepth',mindepth,'reshapex',10);
xlabel('Salinity')
ylabel('Potential temperature')

if isnan(xmin+xmax+ymin+ymax)==0
    set(gca,'xlim',[xmin,xmax],'ylim',[ymin,ymax])
end

subplot(4,8,[4:8,12:16,20:24,28:32])
hold on, grid on, box on
mindepth=PARAM.DEPTH_ZOOM;
[xmin,xmax,ymin,ymax]=plot_prof_hist(CTD,F,'tpot','psal',index_closest_dist,index_closest_date,PARAM,'mindepth',mindepth,'reshapex',20,'reshapey',20);
plot_prof_closest(CTD,F,'tpot','psal',index_closest_dist,index_closest_date,PARAM,'reshapex',20,'reshapey',20);


if isnan(xmin+xmax+ymin+ymax)==0
    set(gca,'xlim',[xmin,xmax],'ylim',[ymin,ymax])
end

%set(gca,'ylim',[2 8])
xlabel('Salinity')
ylabel('Potential temperature')
theax=gca;

if strcmp(PARAM.DATATYPE,'raw') == 1
    hs=suptitle({[num2str(wmonum) ' - Cycle ' num2str(profnum) ' - Raw - Date Argo profile ' datestr(F.dates.data,'dd-mmm-yyyy') ];...
        ['Dates historicals profiles ' datestr(CTD.dates.data(index_closest_dist)) ' (magenta) and ' datestr(CTD.dates.data(index_closest_date)) ' (blue)']});
elseif  strcmp(PARAM.DATATYPE,'adj') == 1
    hs=suptitle({[num2str(wmonum) ' - Cycle ' num2str(profnum) ' - Adj - Date Argo profile ' datestr(F.dates.data,'dd-mmm-yyyy')];...
        ['Dates historicals profiles ' datestr(CTD.dates.data(index_closest_dist)) ' (magenta) and ' datestr(CTD.dates.data(index_closest_date)) ' (blue)']});
end
hs.FontWeight='bold';

%set(hs,'fontsize',14);
colorbar(theax)
% thedate_CTD=datevec(CTD.dates.data);
% thedate_F=datevec(F.dates.data);
% mindate = min(thedate_CTD(:,1));
% mindate = min(mindate,thedate_F(:,1));
% maxdate = max(thedate_CTD(:,1));
% maxdate = max(maxdate,thedate_F(:,1));
% theyear_str=cellstr(num2str([mindate:maxdate]'));
% itick=get(hc,'YTick')
% set(hc,'YTickLabel',theyear_str(itick));
%zone_visu=[minlat-10 maxlat+10 minlon-10 maxlon+10];
subplot(3,8,[17:18])
plot_hist_loc(zone_visu,reso,proj,CTD,F,index_closest_date);
msg='ok';

if PARAM.CHECK_REF==1
    fclose(fw);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_hist_loc(zone_visu,reso,proj,CTD,F,index_closest_date)


fct_pltmap_bathy(zone_visu,reso,proj);
hold on

h0=m_plot(CTD.longitude.data,CTD.latitude.data,'k.');set(h0,'markersize',8);
h0 = m_plot(CTD.longitude.data(index_closest_date),CTD.latitude.data(index_closest_date),'b*'); set(h0,'markersize',12);
h1 = m_plot(F.longitude.data,F.latitude.data,'ro'); set(h1,'markersize',10,'markerfacecolor','r');
h2 = m_plot(CTD.longitude.data(1),CTD.latitude.data(1),'m*'); set(h2,'markersize',12); %closest in space

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_prof_with_flag(F,param1,param2)

tabcol=[0 1 0;1 0.5 0.5;1 0.5 0;1 0 0];
test=check_isfillval_prof(F,param2);
if test.(param2)==0
    hqp=plot(F.(param2).data,F.(param1).data,'g+-','LineWidth',2);
    
    for vqc=2:4
        iqc=find(F.([param2 '_qc']).data == vqc);
        hqc=plot(F.(param2).data(iqc),F.(param1).data(iqc),'+','MarkerSize',3);
        set(hqc,'color',tabcol(vqc,:));
    end
end
return
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_prof_closest(CTD,F,param1,param2,index_closest_dist,index_closest_date,PARAM,varargin);

n=length(varargin);

if n/2~=floor(n/2)
    error('check the imput arguments')
end

f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);

%default
resahpex=1;
reshapey=1;
mindepth=0;

if isfield(s,'reshapex')==1;reshapex=s.reshapex;end;
if isfield(s,'reshapey')==1;reshapey=s.reshapey;end;
if isfield(s,'mindepth')==1;mindepth=s.mindepth;end;


plot(CTD.(param2).data(index_closest_dist,:),CTD.(param1).data(index_closest_dist,:),'m','LineWidth',2);
plot(CTD.(param2).data(index_closest_date,:),CTD.(param1).data(index_closest_date,:),'b','LineWidth',2);

plot_prof_with_flag(F,param1,param2)

if strcmp(param1,'pres')
    set(gca,'ydir','reverse')
end


% xmin = floor(min(min(CTD.(param2).data(:),min(F.(param2).data(:))))*reshapex)/reshapex;
% xmax = ceil(max(max(CTD.(param2).data(:),max(F.(param2).data(:))))*reshapex)/reshapex;
ideC = find(CTD.pres.data>=mindepth);
ideF = find(F.pres.data>mindepth);


if isempty (ideC)==0&isempty(ideF)==0
    xmin = floor(min(min(CTD.(param2).data(ideC),min(F.(param2).data(ideF))))*reshapex)/reshapex;
    xmax = ceil(max(max(CTD.(param2).data(ideC),max(F.(param2).data(ideF))))*reshapex)/reshapex;
 
if ~isnan(xmin)&~isnan(xmax)
    set(gca,'xlim',[xmin xmax]);
end
    %set(gca,'xlim',[xmin xmax]);
    %if ~strcmp(param1,'pres')
    ymin = floor(min(min(F.(param1).data(ideF)))*reshapey)/reshapey;
    ymax = ceil(max(max(F.(param1).data(ideF)))*reshapey)/reshapey;
    if ~isnan(ymin)&~isnan(ymax)
    set(gca,'ylim',[ymin ymax]);
    end

    %set(gca,'ylim',[ymin ymax]);
elseif isempty (ideC)==0
    xmin = floor(min(CTD.(param2).data(ideC))*reshapex)/reshapex;
    xmax = ceil(max(CTD.(param2).data(ideC))*reshapex)/reshapex;
    
if ~isnan(xmin)&~isnan(xmax)
    set(gca,'xlim',[xmin xmax]);
end
    %set(gca,'xlim',[xmin xmax]);
    %if ~strcmp(param1,'pres')
    
    ymin = floor(min(CTD.(param1).data(ideC))*reshapey)/reshapey;
    ymax = ceil(max(CTD.(param1).data(ideC))*reshapey)/reshapey;
    if ~isnan(ymin)&~isnan(ymax)
    set(gca,'ylim',[ymin ymax]);
    end

    %set(gca,'ylim',[ymin ymax]);
else
    xmin=NaN;xmax=NaN;ymin=NaN;ymax=NaN;
end
return

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xmin,xmax,ymin,ymax]=plot_prof_hist(CTD,F,param1,param2,index_closest_dist,index_closest_date,PARAM,varargin);

n=length(varargin);

if n/2~=floor(n/2)
    error('check the imput arguments')
end

f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);

%default
resahpex=1;
reshapey=1;
mindepth=0;

if isfield(s,'reshapex')==1;reshapex=s.reshapex;end;
if isfield(s,'reshapey')==1;reshapey=s.reshapey;end;
if isfield(s,'mindepth')==1;mindepth=s.mindepth;end;


thedate_CTD=datevec(CTD.dates.data);
thedate_F=datevec(F.dates.data);

mindate = min(thedate_CTD(:,1));
mindate = min(mindate,thedate_F(:,1));
maxdate = max(thedate_CTD(:,1));
maxdate = max(maxdate,thedate_F(:,1));
%c=cmocean('thermal',maxdate-mindate+1);
l_color=maxdate-mindate+1;
d_color=floor(l_color/4);
%c1=gray(l_color*2);
c1=cmocean('dense',l_color*2);
c=c1(d_color:d_color+l_color,:);


%c=c1([floor(length(c1)/4):3*ceil(length(c1)/4)+4],:)
colormap(c)
theyear=[mindate:maxdate];


for i=1:length(CTD.dates.data)
    idate=find(theyear==thedate_CTD(i,1));
    %thecolomark=c(idate,:);
    %plot(CTD.(param2).data(i,:),CTD.(param1).data(i,:),'.k','Markersize',10);
    if strcmp(PARAM.REFERENCE,'argo') == 1
        scatter(CTD.(param2).data(i,:),CTD.(param1).data(i,:),7,repmat(thedate_CTD(i,1),[1,size(CTD.(param2).data(i,:),2)]),'*')
    elseif strcmp(PARAM.REFERENCE,'ctd') == 1
        scatter(CTD.(param2).data(i,:),CTD.(param1).data(i,:),3,repmat(thedate_CTD(i,1),[1,size(CTD.(param2).data(i,:),2)]),'*')
    end
    hold on
    %h=plot(CTD.(param2).data(i,:),CTD.(param1).data(i,:),'.k','Markersize',10);
    %set(h,'color',c(idate,:));
end


plot_prof_with_flag(F,param1,param2)

if strcmp(param1,'pres')
    set(gca,'ydir','reverse')
end

ideC = find(CTD.pres.data>=mindepth);
ideF = find(F.pres.data>mindepth);
if isempty (ideC)==0&isempty(ideF)==0
    xmin = floor(min(min(CTD.(param2).data(ideC),min(F.(param2).data(ideF))))*reshapex)/reshapex;
    xmax = ceil(max(max(CTD.(param2).data(ideC),max(F.(param2).data(ideF))))*reshapex)/reshapex;
   
if ~isnan(xmin)&~isnan(xmax)
    set(gca,'xlim',[xmin xmax]);
end
    %set(gca,'xlim',[xmin xmax]);
    mindepth;
    if (mindepth)==0
        ymin = floor(min(min(CTD.(param1).data(ideC),min(F.(param1).data(ideF))))*reshapey)/reshapey;
        ymax = ceil(max(max(CTD.(param1).data(ideC),max(F.(param1).data(ideF))))*reshapey)/reshapey;
    else
        ymin = floor(min(min(F.(param1).data(ideF)))*reshapey)/reshapey;
        ymax = ceil(max(max(F.(param1).data(ideF)))*reshapey)/reshapey;
    end
    %ymin = floor(min(min(CTD.(param1).data(ideC),min(F.(param1).data(ideF))))*reshapey)/reshapey;
    %ymax = ceil(max(max(CTD.(param1).data(ideC),max(F.(param1).data(ideF))))*reshapey)/reshapey;
    [ymin ymax];
    if ~isnan(ymin)&~isnan(ymax)
    set(gca,'ylim',[ymin ymax]);
    end

   % set(gca,'ylim',[ymin ymax]);
else
    xmin=NaN;xmax=NaN;ymin=NaN;ymax=NaN;
end

return

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PV=compute_PV(LONG,LAT)

la_grid_long1=LONG;
gg=find(LONG>180);
la_grid_long1(gg)=LONG(gg)-360; % m_tbase inputs longitudes from 0 to +/- 180
xlim=  [min(la_grid_long1)-1, max(la_grid_long1)+1];
if xlim(1)<-180; xlim(1)=-180;end;
if xlim(2)>180; xlim(2)=180;end;
ylim=  [min(LAT)-1, max(LAT)+1];
if ylim(1)<-90; ylim(1)=-90;end;
if ylim(2)>90; ylim(2)=90;end;
m_proj('mercator','long', xlim, 'lat', ylim);
[elev,x,y] = m_tbase( [min(la_grid_long1)-1, max(la_grid_long1)+1, min(LAT)-1, max(LAT)+1] );
la_grid_Z = -interp2( x,y,elev, la_grid_long1, LAT, 'linear'); % -ve bathy values

PV= (2*7.292*10^-5.*sin(LAT.*pi/180))./la_grid_Z;

jj=find(PV==0);
PV(jj)=1*10^-5;
return


