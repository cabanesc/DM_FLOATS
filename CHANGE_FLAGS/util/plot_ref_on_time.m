function plot_ref_on_time(DATA_BASE,CONFIG,Param1,Param,varargin)
% -========================================================
%   USAGE : plot_ref_on_time(DATA_BASE,CONFIG,Param1,Param,varargin)
%   PURPOSE : trace les salinité de reference  a un niveau theta donne, en fonction du temps
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
accepted_Param1={'','psal_mean','thedates','pres_mean','longitude','latitude','tpot_mean'};
accepted_Param={'','psal_mean','thedates','pres_mean','longitude','latitude','tpot_mean'};
if ismember(Param,accepted_Param)==0
    disp(['Accepted param: ' accepted_Param])
    error(['Unreconized Param ' Param])
end
thefield = fieldnames(s);
for k=1:length(thefield)
    if ismember(thefield{k},accepted_options)==0
        disp(['Accepted options: ' accepted_options])
        error(['Unreconized option ' thefield{k}])
    end
end
if isfield(s,'MarkerSize')==0
    s.MarkerSize=20; %default
end
if isfield(s,'MarkerSpec')==0
    s.MarkerSpec='o'; %default
end
if isempty(Param)==1&~isfield(s,'MarkerFaceColor')
    s.MarkerFaceColor='m'; % default
end

ModeClim=get(gca,'CLimMode')
ht=get(gca,'Title');
thetitle=get(ht,'string')

thestr='';

%keyboard
if isempty(strfind(Param1,'psal_mean'))==0|isempty(strfind(Param,'psal_mean'))==0
   thestr=[thestr 'Mean salinity at'];
elseif isempty(strfind(Param1,'tpot_mean'))==0|isempty(strfind(Param,'tpot_mean'))==0
   thestr=[thestr 'Mean potential temperature at'];
else
  thestr=[thestr Param ' at'];
end
if CONFIG.OnTheta==1
   thestr=[thestr ' theta: ' num2str(CONFIG.TPOT_MIN) '-' num2str(CONFIG.TPOT_MAX)];
else
   thestr=[thestr ' P: ' num2str(CONFIG.P_MIN) '-' num2str(CONFIG.P_MAX)];
end

if ~isempty(thetitle)
        if isempty(strfind(thetitle,DATA_BASE))&(isempty(strfind(thestr,'salinity'))&isempty(strfind(thetitle,'psal_mean'))|isempty(strfind(thestr,'temperature'))&isempty(strfind(thetitle,'tpot_mean')))
            title([thestr ' from ARGO and CTD reference databases']);
        else
            cla(gca)
            ModeClim='auto';
            title([thestr ' from ' DATA_BASE ' reference database']);
        end
else
        title([thestr ' from ' DATA_BASE ' reference database']);
end
    


h=get(gca,'Children');
thetexth=findobj(h,'type','text','rotation',270);
thetext=get(thetexth,'string')

if ~isempty(thetext)
    if isempty(strfind(thetext,upper(Param)))
        cla(gca)
        ModeClim='auto';
    else
        %keyboard
        set(thetexth,'string','')
    end
end
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
    if CONFIG.OnTheta==1
        [PLOT,CTD] = find_psal_on_theta(CTD, CONFIG.TPOT_MIN, CONFIG.TPOT_MAX, CONFIG.MIN_DEPTH,PLOT);
        %Param1='psal_mean';
    else
        [PLOT,CTD] = find_psal_on_z(CTD, CONFIG.P_MIN, CONFIG.P_MAX,PLOT);
        %Param1='tpot_mean';
    end
    CTD.thedates.data=((CTD.dates.data-datenum('19700101','yyyymmdd'))/365.25+1970);
    ip=find(CTD.latitude.data>=CONFIG.VEC_REG(3)&CTD.latitude.data<=CONFIG.VEC_REG(4)&CTD.longitude.data<=CONFIG.VEC_REG(2)&CTD.longitude.data>=CONFIG.VEC_REG(1)&CTD.thedates.data>=CONFIG.YEARDEB&CTD.thedates.data<CONFIG.YEARFIN+1);
    if isempty(Param)==0
        if ismember(s.MarkerSpec,{'.','*','x','+'})
            scatter(CTD.thedates.data(ip),CTD.(Param1).data(ip),s.MarkerSize,CTD.(lower(Param)).data(ip),s.MarkerSpec)
        else
            scatter(CTD.thedates.data(ip),CTD.(Param1).data(ip),s.MarkerSize,CTD.(lower(Param)).data(ip),s.MarkerSpec,'filled')
        end
    else
        if ismember(s.MarkerSpec,{'.','*','x','+'})
            scatter(CTD.thedates.data(ip),CTD.(Param1).data(ip),s.MarkerSize,[s.MarkerFaceColor s.MarkerSpec])
        else
            scatter(CTD.thedates.data(ip),CTD.(Param1).data(ip),s.MarkerSize,[s.MarkerFaceColor s.MarkerSpec],'filled')
        end
    end
    if isempty(ip)==0
        themin=min(themin,min(CTD.(lower(Param)).data(ip)));
        themax=max(themax,max(CTD.(lower(Param)).data(ip)));
    end
    ylabel(Param1,'interpreter','none')
end
%keyboard
caxis([themin,themax])
colorbar
b=get(gca,'YLim');
%axis([CONFIG.YEARDEB CONFIG.YEARFIN b(1) b(2)])
a=get(gca,'XLim');

text(a(2)+(a(2)-a(1))/4.7,b(1)+2/3*(b(2)-b(1)),upper(Param),'Rotation',-90,'FontSize',11,'interpreter','none')
