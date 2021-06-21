
function plot_diagnostics_ow_figure9_2( pn_float_dir, pn_float_name, po_system_configuration,dacname )

%
% Annie Wong, 14 June 2011
% Breck Owens, October 2006
%--------------------------------------------------------------------------

%pn_float_dir='uw/';
%pn_float_name='R5902134';
%po_system_configuration = load_configuration( 'ow_config.txt' );


%close all

% modify pn_float_name for title if the name contains '_' ------------

ii=find(pn_float_name=='_');
if(isempty(ii)==0)
  title_floatname = strcat( pn_float_name(1:ii-1), '\_', pn_float_name(ii+1:length(pn_float_name)));
else
  title_floatname = pn_float_name;
end


% load data from /float_source, /float_mapped, /float_calib, and config --------------

lo_float_source_data = load(fullfile( po_system_configuration.FLOAT_SOURCE_DIRECTORY, pn_float_dir, ...
  strcat( pn_float_name, po_system_configuration.FLOAT_SOURCE_POSTFIX ) ) );

Co.longitude.data=lo_float_source_data.LONG;
Co.latitude.data=lo_float_source_data.LAT;


Co=shiftEW(Co,'longitude','grwch');

PROFILE_NO = lo_float_source_data.PROFILE_NO;
LAT  = Co.latitude.data;
LONG = Co.longitude.data;
PRES = lo_float_source_data.PRES;
TEMP = lo_float_source_data.TEMP;
PTMP = lo_float_source_data.PTMP;
SAL  = lo_float_source_data.SAL;

if(isempty(find(isnan(PRES)==0))==0) % if no data exists, terminate here, no plots will be produced

lo_float_mapped_data = load( fullfile( po_system_configuration.FLOAT_MAPPED_DIRECTORY, pn_float_dir, ...
  strcat( po_system_configuration.FLOAT_MAPPED_PREFIX, pn_float_name, po_system_configuration.FLOAT_MAPPED_POSTFIX ) ) ) ;

mapped_sal  = lo_float_mapped_data.la_mapped_sal;
mapsalerrors = lo_float_mapped_data.la_mapsalerrors;
la_ptmp = lo_float_mapped_data.la_ptmp;
selected_hist =lo_float_mapped_data.selected_hist;

if(isempty(find(isnan(mapped_sal)==0))==0) % if mapping exists, terminate here, no plots will be produced

lo_float_calib_data = load( fullfile( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, ...
  strcat( po_system_configuration.FLOAT_CALIB_PREFIX, pn_float_name, po_system_configuration.FLOAT_CALIB_POSTFIX ) ) );

cal_SAL = lo_float_calib_data.cal_SAL;
cal_SAL_err = lo_float_calib_data.cal_SAL_err;
pcond_factor = lo_float_calib_data.pcond_factor;
pcond_factor_err = lo_float_calib_data.pcond_factor_err;

lo_float_calseries = load( fullfile( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir, ...
  strcat( po_system_configuration.FLOAT_CALSERIES_PREFIX , pn_float_name, po_system_configuration.FLOAT_MAPPED_POSTFIX ) ) );

use_theta_gt = lo_float_calseries.use_theta_gt;
use_theta_lt = lo_float_calseries.use_theta_lt;
use_pres_gt = lo_float_calseries.use_pres_gt;
use_pres_lt = lo_float_calseries.use_pres_lt;
use_percent_gt = lo_float_calseries.use_percent_gt;

% load the station by station fits
load(fullfile( po_system_configuration.FLOAT_CALIB_DIRECTORY, pn_float_dir,...
   strcat( po_system_configuration.FLOAT_CALIB_PREFIX, pn_float_name, po_system_configuration.FLOAT_CALIB_POSTFIX ) ),'-regexp','^sta')


% plot the float locations (figure 1) -----------------------

load( fullfile( po_system_configuration.CONFIG_DIRECTORY, po_system_configuration.CONFIG_COASTLINES ), 'coastdata_x', 'coastdata_y' );


%keyboard
[m,n] = size(PRES);

%figure
subplot(2,2,2)


Topo_ficin='/home/lpoargo1/DMARGO/OW/TOPO/topo.onetenthdeg.atl.nc';
Topo=read_netcdf_allthefile(Topo_ficin);
hold on

%load colormap.mat

contour(Topo.lon.data,Topo.lat.data,Topo.topo.data,[-1000 -1000],'LineColor',[0.6 0.5 0.4]);
contour(Topo.lon.data,Topo.lat.data,Topo.topo.data,[-2000 -2000],'LineColor',[0.6 0.5 0.4]);
contour(Topo.lon.data,Topo.lat.data,Topo.topo.data,[-2500 -2500],'LineColor',[0.8 0.7 0.6]);

contour(Topo.lon.data,Topo.lat.data,Topo.topo.data,[0 0],'LineColor',[0 0 0],'LineWidth',1.5);        
grid on
%colormap(w)
box on

%  set(gcf,'defaultaxeslinewidth',2)
%  set(gcf,'defaultlinelinewidth',2)
%  set(gcf,'defaultaxesfontsize',12)

colormap(jet(n));
c=colormap;

x=[];
y=[];
if(isempty(selected_hist)==0)
    HIST.longitude.data=selected_hist(:,1);
    HIST.latitude.data=selected_hist(:,2);
    HIST=shiftEW(HIST,'longitude','grwch');
    x=HIST.longitude.data;
    y=HIST.latitude.data;
end



plot(LONG,LAT,'r-','LineWidth',1);
hold on
plot(x,y,'.','Color',[0.5 0.5 0.5],'MarkerSize',8)

%legend('float','historical points','Location','Best')
%plot(coastdata_x,coastdata_y,'k.-');
for i=1:n
  h=plot(LONG(i),LAT(i),'+','MarkerSize',4);
  set(h,'color',c(i,:));
 % j=text(LONG(i),LAT(i),int2str(PROFILE_NO(i)));
  %set(j,'color',c(i,:),'fontsize',12,'hor','cen');
end

minx = max(min(min(x),min(LONG))-5,-100);
maxx = max(max(x),max(LONG))+5;
miny = min(min(y),min(LAT))-5;
maxy=  max(max(y),max(LAT))+5;

if maxy<75; maxy=maxy;end

if  maxy>=75; maxy=70;end


axis([minx maxx miny maxy])

set(gca,'FontSize',11)
xlabel('Longitude');
ylabel('Latitude');
title( strcat(' CTD REF : ', po_system_configuration.CONFIGURATION_FILE(1:end-4) ),'interpreter','none' );




% calibration curve (figure 3) --------------------------

if(isempty(find(isnan(cal_SAL)==0))==0) % if no cal exists, terminate here, no plot will be produced

Soffset=cal_SAL-SAL;
avg_Soffset=NaN.*ones(1,n);
avg_Soffset_err=NaN.*ones(1,n);
Staoffset=sta_SAL-SAL;
avg_Staoffset=NaN.*ones(1,n);
avg_Staoffset_err=NaN.*ones(1,n);

for i=1:n
   ii=[];
   ii=find(isnan(Soffset(:,i))==0);
   if ~isempty(ii)
       avg_Soffset(i)=mean(Soffset(ii,i));
       avg_Soffset_err(i)=mean(cal_SAL_err(ii,i));
   else
       avg_Soffset(i) = NaN;
       avg_Soffset_err(i) = NaN;
   end
   ii=find(isnan(Staoffset(:,i))==0);
   if ~isempty(ii)
       avg_Staoffset(i)=mean(Staoffset(ii,i));
       avg_Staoffset_err(i)=mean(sta_SAL_err(ii,i));
   else
       avg_Staoffset(i) = NaN;
       avg_Staoffset_err(i) = NaN;
   end
end
 


subplot(2,1,2)
%  plot(PROFILE_NO, avg_Soffset, 'b-');
%  hold on
%  plot(PROFILE_NO, avg_Soffset, 'g-');
% Plot station by station fit
ok = find(isfinite(avg_Staoffset));
plot(PROFILE_NO(ok), avg_Staoffset(ok), 'r+','MarkerSize',10);

%  legend('2 x cal error','1 x cal error','1-1 profile fit', 'DM correction', 'Location', 'SouthOutside','orientation','horizontal');
%  errorbar(PROFILE_NO, avg_Soffset, 2*avg_Soffset_err,'b')
%  errorbar(PROFILE_NO, avg_Soffset, avg_Soffset_err,'go-')
%  plot(PROFILE_NO(ok), avg_Staoffset(ok), 'r-');

%axis([ 0, max(PROFILE_NO)+1, min([avg_Soffset-avg_Soffset_err,0])-.02, max([avg_Soffset+avg_Soffset_err,0])+.02 ])
%plot( [0, max(PROFILE_NO)+1], [0,0], 'k-')
axis tight
set(gca,'FontSize',10)
%  xlabel(['float profile number (date of the 1 cycle: ' thedate_1cycle ]);
%  ylabel('\Delta S')
%  title( strcat(' vertically-averaged salinity (PSS-78) additive correction \Delta S (red) with errors and DM correction (magenta)') );
%  grid on

axis tight
v=axis;
axis([v(1) v(2) min(v(3),-0.01) max(v(4),0.01)])

v=axis;
%keyboard
hold on
for ip=1:length(PROFILE_NO)
h=plot( PROFILE_NO(ip), v(3)-0.001,'s','MarkerSize',4);
  %set(h,'color',c(ip,:));
  set(h,'color',c(ip,:),'MarkerFaceColor',c(ip,:))
end

end %if(isempty(find(isnan(cal_SAL)==0))==0) ---------------


drawnow
set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.75,8,9.5]);
print('-depsc ', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_bilan.eps'));


disp(' ')
disp(' ')



else %if(isempty(find(isnan(mapped_sal)==0))==0) -----------
drawnow
set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.75,8,9.5]);
print('-depsc ', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_bilan.eps'));
end

end %if(isempty(find(isnan(PRES)==0))==0) ------------------


