function plot_ref_psal(DATA_BASE,CONFIG,varargin)
% -========================================================
%   USAGE : plot_ref_psal(DATA_BASE,CONFIG,varargin)
%   PURPOSE : trace les diagramme theta/S des donnees de reference
% -----------------------------------
%   INPUT :
%     IN1   (class)  -comments-
%             additional description
%     IN2   (class)  -comments-
%
%   OPTIONNAL INPUT :
%    OPTION1  (class)  -comments-
% -----------------------------------
%   OUTPUT :
%     OUT1   (class)  -comments-
%             additional description
%     OUT2   (class)  -comments-
%             additional description
% -----------------------------------
%   HISTORY  : created (2009) ccabanes
%            : modified (yyyy) byxxx 
%   CALLED SUBROUTINES: none
% ========================================================
n=length(varargin);

if n/2~=floor(n/2)
    error('check the imput arguments')
end

f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);
accepted_options={'MarkerSize','MarkerFaceColor','MarkerSpec'};
accepted_Param={'','psal_mean','thedates','pres_mean','longitude','latitude'};
%  if ismember(Param,accepted_Param)==0
%     error(['Unreconized Param ' Param])
%  end
thefield = fieldnames(s);
for k=1:length(thefield)
    if ismember(thefield{k},accepted_options)==0
    error(['Unreconized option ' thefield{k}])
    end
end
if isfield(s,'MarkerSize')==0
s.MarkerSize=20; %default
end
if isfield(s,'MarkerSpec')==0
s.MarkerSpec='o'; %default
end
if ~isfield(s,'MarkerFaceColor')
   s.MarkerFaceColor='m'; % default
end

ModeClim=get(gca,'CLimMode')
ht=get(gca,'Title');
thetitle=get(ht,'string')



%keyboard
if ~isempty(thetitle)
 if isempty(strfind(thetitle,DATA_BASE))
        title(['Salinity from ARGO and CTD reference databases']);
 else
        %cla(gca)
        ModeClim='auto';
        title(['Theta/S diagram from ' DATA_BASE ' reference database']);
 end
else
title(['Salinity from ' DATA_BASE ' reference database']);
end

hold on
box on
grid on

region.lonmin=CONFIG.VEC_REG(1);
region.lonmax=CONFIG.VEC_REG(2);
region.latmin=CONFIG.VEC_REG(3);
region.latmax=CONFIG.VEC_REG(4);
DATA_BASE
% Identification des boites wmo à vérifier
switch DATA_BASE
    case 'CTD'
        load([CONFIG.pathwmobox 'constants/wmo_boxes_ctd.mat']);
        wmo_square = find_wmoboxes( region, la_wmo_boxes);
        la_wmo_num=la_wmo_boxes(ismember(la_wmo_boxes(:,1),wmo_square)&la_wmo_boxes(:,2)==1,:);
        rep='climatology/historical_ctd/ctd_';
    case 'ARGO'
        load([CONFIG.pathwmobox 'constants/wmo_boxes_argo.mat']);
        wmo_square = find_wmoboxes( region, la_wmo_boxes);
        la_wmo_num=la_wmo_boxes(ismember(la_wmo_boxes(:,1),wmo_square)&la_wmo_boxes(:,4)==1,:);
        rep='climatology/argo_profiles/argo_';
end

[nbox,null]=size(la_wmo_num);

if strcmp(CONFIG.VEC_PSAL_AUTO,'on')==0
caxis(CONFIG.VEC_PSAL)
end


if isfield(CONFIG,'YEARDEB')==0
CONFIG.YEARDEB=-9999;
end
if isfield(CONFIG,'YEARFIN')==0
CONFIG.YEARFIN=9999;
end
alldates=[];
switch ModeClim
case 'auto'
themin=999999;
themax=-999999;
case 'manual'
color_lim=get(gca,'Clim');
themin=color_lim(1);
themax=color_lim(2);
end

hold on
box on
grid on
for ibox=1:nbox % boucle sur toutes les boites
    clear CTD
    PLOT=[];
    wmo_file= [CONFIG.pathwmobox rep num2str(la_wmo_num(ibox,1))];
    [CTD] = read_wmo(wmo_file,CONFIG.ctd_to_exclude,DATA_BASE);
     CTD.tpot_min.data=min(CTD.tpot.data');
     CTD.pres_max.data=max(CTD.pres.data');
    [PLOT,CTD] = find_psal_on_theta(CTD, CONFIG.TPOT_MIN, CONFIG.TPOT_MAX, CONFIG.MIN_DEPTH,PLOT);
    CTD.thedates.data=((CTD.dates.data-datenum('19700101','yyyymmdd'))/365.25+1970);
    ip=find(CTD.latitude.data>=CONFIG.VEC_REG(3)&CTD.latitude.data<=CONFIG.VEC_REG(4)&CTD.longitude.data<=CONFIG.VEC_REG(2)&CTD.longitude.data>=CONFIG.VEC_REG(1)&CTD.thedates.data>=CONFIG.YEARDEB&CTD.thedates.data<CONFIG.YEARFIN);
    
    for k=1:length(ip)
        
        plot(CTD.psal.data(ip(k),:)',CTD.pres.data(ip(k),:)','Color',s.MarkerFaceColor);
       
    end
%      if isempty(ip)==0
%      themin=min(themin,min(CTD.(lower(Param)).data(ip)));
%      themax=max(themax,max(CTD.(lower(Param)).data(ip)));
%      end
end 
set(gca,'Ydir','reverse')
%keyboard
%  caxis([themin,themax])
%  colorbar
%  a=get(gca,'XLim');
%  b=get(gca,'YLim');
%  text(a(2)+(a(2)-a(1))/4.7,b(1)+2/3*(b(2)-b(1)),upper(Param),'Rotation',-90,'FontSize',11)
