function plot_ref_on_map(DATA_BASE,topo,CONFIG,Param,varargin)
% -========================================================
%   USAGE : plot_ref_on_map(DATA_BASE,topo,CONFIG,Param,varargin)
%   PURPOSE : trace les données de reference sur une carte a un niveau theta donne ou P donne
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
accepted_Param={'','psal_mean','thedates','pres_mean','tpot_min','pres_max','tpot_mean'};
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
%keyboard
if ~isempty(thetitle)
    switch Param
        case 'psal_mean'
            if  ~isempty(strfind(thetitle,'salinity'))
                if isempty(strfind(thetitle,DATA_BASE))
                    disp('ok')
                    if CONFIG.ONTHETA==1
                        title(['Mean salinity at theta:' num2str(CONFIG.TPOT_MIN) '-' num2str(CONFIG.TPOT_MAX) ' from ARGO and CTD reference databases']);
                    else
                        title(['Mean salinity at P:' num2str(CONFIG.P_MIN) '-' num2str(CONFIG.P_MAX) ' from ARGO and CTD reference databases']);
                    end
                end
            else
                cla(gca)
                ModeClim='auto';thetitle='';
            end
         case 'tpot_mean'
            if  ~isempty(strfind(thetitle,'temperature'))
                if isempty(strfind(thetitle,DATA_BASE))
                    disp('ok')
                    if CONFIG.ONTHETA==1
                        title(['Mean potential temperature at theta:' num2str(CONFIG.TPOT_MIN) '-' num2str(CONFIG.TPOT_MAX) ' from ARGO and CTD reference databases']);
                    else
                        title(['Mean potential temperature at P:' num2str(CONFIG.P_MIN) '-' num2str(CONFIG.P_MAX) ' from ARGO and CTD reference databases']);
                    end
                end
            else
                cla(gca)
                ModeClim='auto';thetitle='';
            end   
        case 'pres_mean'
            if  ~isempty(strfind(thetitle,'pressure'))
                if isempty(strfind(thetitle,DATA_BASE))
                    title(['Mean pressure at theta:' num2str(CONFIG.TPOT_MIN) '-' num2str(CONFIG.TPOT_MAX) ' from ARGO and CTD reference databases']);
                end
            else
                cla(gca)
                ModeClim='auto';thetitle='';
            end
        case 'thedates'
            if  ~isempty(strfind(thetitle,'Dates'))
                if isempty(strfind(thetitle,DATA_BASE))
                    title(['Dates of measurements,  ARGO and CTD reference databases']);
                end
            else
                cla(gca)
                ModeClim='auto';thetitle='';
            end
        case 'tpot_min'
            if  ~isempty(strfind(thetitle,'Minimum'))
                if isempty(strfind(thetitle,DATA_BASE))
                    title(['Minimum potential temperature on the profiles (ARGO and CTD reference database)'])
                end
            else
                cla(gca)
                ModeClim='auto';thetitle='';
            end
        case 'pres_max'
            if  ~isempty(strfind(thetitle,'Maximum'))
                if isempty(strfind(thetitle,DATA_BASE))
                    title(['Maximum pressure on the profiles (ARGO and CTD reference database)'])
                end
            else
                cla(gca)
                ModeClim='auto';thetitle='';
            end
    end
end
if isempty(thetitle)
    switch Param
        case 'psal_mean'
            if CONFIG.ONTHETA==1
                title(['Mean salinity at theta:' num2str(CONFIG.TPOT_MIN) '-' num2str(CONFIG.TPOT_MAX) ' from ' DATA_BASE ' reference database']);
            else
                title(['Mean salinity at P:' num2str(CONFIG.P_MIN) '-' num2str(CONFIG.P_MAX) ' from ' DATA_BASE ' reference database']);
            end
        case 'pres_mean'
            title(['Mean pressure at theta:' num2str(CONFIG.TPOT_MIN) '-' num2str(CONFIG.TPOT_MAX) ' from ' DATA_BASE ' reference database']);
        case 'thedates'
            title(['Dates of measurements, ' DATA_BASE ' reference database']);
        case 'tpot_min'
            title(['Minimum potential temperature on the profiles (' DATA_BASE ' reference database)'])
        case 'pres_max'
            title(['Maximum pressure on the profiles (' DATA_BASE ' reference database)'])
    end
end
%keyboard

hold on
box on
grid on

contour_bathy(topo,CONFIG.VEC_REG)

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


if isfield(CONFIG,'YEAR_MIN')==0
    CONFIG.YEAR_MIN=-9999;
end
if isfield(CONFIG,'YEAR_MAX')==0
    CONFIG.YEAR_MAX=9999;
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
allsource={''};
for ibox=1:nbox % boucle sur toutes les boites
    clear CTD
    PLOT=[];
    wmo_file= [CONFIG.pathwmobox rep num2str(la_wmo_num(ibox,1))];
    [CTD] = read_wmo(wmo_file,CONFIG.ctd_to_exclude,DATA_BASE);
    CTD.tpot_min.data=min(CTD.tpot.data');
    CTD.pres_max.data=max(CTD.pres.data');
    if CONFIG.ONTHETA==1
        [PLOT,CTD] = find_psal_on_theta(CTD, CONFIG.TPOT_MIN, CONFIG.TPOT_MAX, CONFIG.DEPTH_MIN,PLOT);
    else
        [PLOT,CTD] = find_psal_on_z(CTD, CONFIG.P_MIN, CONFIG.P_MAX,PLOT);
    end
    CTD.thedates.data=((CTD.dates.data-datenum('19700101','yyyymmdd'))/365.25+1970);
    ip=find(CTD.latitude.data>=CONFIG.VEC_REG(3)&CTD.latitude.data<=CONFIG.VEC_REG(4)&CTD.longitude.data<=CONFIG.VEC_REG(2)&CTD.longitude.data>=CONFIG.VEC_REG(1)&CTD.thedates.data>=CONFIG.YEAR_MIN&CTD.thedates.data<CONFIG.YEAR_MAX+1);
    allsource=[allsource unique(strtok(CTD.source.data(ip),'_'))];
    if isempty(Param)==0
        if ismember(s.MarkerSpec,{'.','*','x','+'})
            scatter(CTD.longitude.data(ip),CTD.latitude.data(ip),s.MarkerSize,CTD.(lower(Param)).data(ip),s.MarkerSpec)
        else
            scatter(CTD.longitude.data(ip),CTD.latitude.data(ip),s.MarkerSize,CTD.(lower(Param)).data(ip),s.MarkerSpec,'filled')
        end
    else
        if ismember(s.MarkerSpec,{'.','*','x','+'})
            scatter(CTD.longitude.data(ip),CTD.latitude.data(ip),s.MarkerSize,[s.MarkerFaceColor s.MarkerSpec])
        else
            scatter(CTD.longitude.data(ip),CTD.latitude.data(ip),s.MarkerSize,[s.MarkerFaceColor s.MarkerSpec],'filled')
        end
    end
    if isempty(ip)==0
        if isempty(Param)==0
            themin=min(themin,min(CTD.(lower(Param)).data(ip)));
            themax=max(themax,max(CTD.(lower(Param)).data(ip)));
        end
    end
end
if isempty(Param)==0
    caxis([themin,themax])
    colorbar
end

disp('source utilisées:')
unique(allsource)
