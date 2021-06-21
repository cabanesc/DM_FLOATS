function plot_ref_diag_close_Argo(DATA_BASE,CONFIG,floatname,dacname,cycle,topo,varargin)
% -========================================================
%   USAGE : plot_ref_diag(DATA_BASE,CONFIG,varargin)
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
addpath('/home5/pharos/argo/DMARGO/OW/VERSION_1_1_5_geovide/matlab_codes')
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
        title(['Theta/S diagram from ARGO and CTD reference databases']);
    else
        %cla(gca)
        ModeClim='auto';
        title(['Theta/S diagram from ' DATA_BASE ' reference database']);
    end
else
    title(['Theta/S diagram from ' DATA_BASE ' reference database']);
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

filename=[CONFIG.DIR_FTP  dacname '/' floatname '/' floatname '_prof.nc'];
F=read_netcdf_allthefile(filename);
F = replace_fill_bynan(F);
F = format_flags_char2num(F);
F.psal.data(F.psal_qc.data>2)=NaN;
F.tpot.data = sw_ptmp(F.psal.data,F.temp.data,F.pres.data,0);

cy=find(F.cycle_number.data==cycle);
LONG=F.longitude.data(cy);
LAT=F.latitude.data(cy);
thedate = datevec((F.juld.data+datenum('19500101','yyyymmdd')));
siz=size(thedate,1);
ll=[thedate(:,1),ones(siz,2),zeros(siz,3)];
F.thedates.data = thedate(:,1)+etime(thedate,ll)./(3600*24*365.25);
DATES=F.thedates.data(cy);
if( isnan(LONG)==0 & isnan(LAT)==0)
    
    if(LONG>180) % m_tbase inputs longitudes from 0 to +/- 180
        LONG1=LONG-360;
    else
        LONG1=LONG;
    end
    m_proj('mercator','long', [min(LONG1)-1, max(LONG1)+1], 'lat', [min(LAT)-1, max(LAT)+1] );
    [elev,x,y] = m_tbase( [min(LONG1)-1, max(LONG1)+1, min(LAT)-1, max(LAT)+1] );
    Z = -interp2( x,y,elev, LONG1, LAT, 'linear'); % -ve bathy values
end
subplot(1,2,1)
contour_bathy(topo,CONFIG.VEC_REG)

for ibox=1:nbox % boucle sur toutes les boites
    
    clear CTD
    PLOT=[];
    wmo_file= [CONFIG.pathwmobox rep num2str(la_wmo_num(ibox,1))];
    [CTD] = read_wmo(wmo_file,CONFIG.ctd_to_exclude,DATA_BASE);
    CTD.tpot_min.data=min(CTD.tpot.data');
    CTD.pres_max.data=max(CTD.pres.data');
    [PLOT,CTD] = find_psal_on_theta(CTD, CONFIG.TPOT_MIN, CONFIG.TPOT_MAX, CONFIG.MIN_DEPTH,PLOT);
    CTD.thedates.data=((CTD.dates.data-datenum('19700101','yyyymmdd'))/365.25+1970);
    
    %keyboard
    la_grid_long= CTD.longitude.data;
    la_grid_long1= CTD.longitude.data;
    gg=find( CTD.longitude.data>180);
    la_grid_long1(gg)= CTD.longitude.data(gg)-360; % m_tbase inputs longitudes from 0 to +/- 180
    la_grid_lat= CTD.latitude.data;
    m_proj('mercator','long', [min(la_grid_long1)-1, max(la_grid_long1)+1], 'lat', [min(la_grid_lat)-1, max(la_grid_lat)+1] );
    [elev,x,y] = m_tbase( [min(la_grid_long1)-1, max(la_grid_long1)+1, min(la_grid_lat)-1, max(la_grid_lat)+1] );
    la_grid_Z = -interp2( x,y,elev, la_grid_long1, la_grid_lat, 'linear'); % -ve bathy values
    PV_float = (2*7.292*10^-5.*sin(LAT.*pi/180))./Z;
    PV_hist = (2*7.292*10^-5.*sin(la_grid_lat.*pi/180))./la_grid_Z;
    if(PV_float==0)PV_float=1*10^-5;end
    jj=find(PV_hist==0);
    PV_hist(jj)=1*10^-5;
    ellipse_large = sqrt( (la_grid_long-LONG).^2./(CONFIG.longitude_large).^2 + (la_grid_lat-LAT).^2./(CONFIG.latitude_large).^2 +...
        ((PV_float-PV_hist)./sqrt( PV_float.^2+PV_hist.^2 )./CONFIG.phi_large).^2 ) ;
    
    ellipse_small = sqrt( (la_grid_long-LONG).^2./(CONFIG.longitude_small).^2 + (la_grid_lat-LAT).^2./(CONFIG.latitude_small).^2 +...
        ((PV_float-PV_hist)./sqrt( PV_float.^2+PV_hist.^2 )./CONFIG.phi_small).^2 ) ;
    ip=find(ellipse_large<1);
    subplot(1,2,1)
    hold on
    box on
    grid on
    scatter(CTD.longitude.data(ip,:),CTD.latitude.data(ip,:),30,s.MarkerFaceColor,'filled')
    scatter(LONG,LAT,30,'m','filled')
    for k=1:length(ip)
        % keyboard
        subplot(1,2,2)
        hold on
        box on
        grid on
        scatter(CTD.psal.data(ip(k),:)',CTD.tpot.data(ip(k),:)',10,repmat(CTD.thedates.data(ip(k)),[1,size(CTD.psal.data(ip(k),:),2)]),'o','filled')
        %plot(CTD.psal.data(ip(k),:)',CTD.tpot.data(ip(k),:)','Color',s.MarkerFaceColor);
        
    end
    %      if isempty(ip)==0
    %      themin=min(themin,min(CTD.(lower(Param)).data(ip)));
    %      themax=max(themax,max(CTD.(lower(Param)).data(ip)));
    %      end
end
subplot(1,2,2)
plot(F.psal.data(cy,:)',F.tpot.data(cy,:)','m','LineWidth',2);

%keyboard
%  caxis([themin,themax])
%  colorbar
%  a=get(gca,'XLim');
%  b=get(gca,'YLim');
%  text(a(2)+(a(2)-a(1))/4.7,b(1)+2/3*(b(2)-b(1)),upper(Param),'Rotation',-90,'FontSize',11)
