
function plot_diagnostics_ow_figure9_3( pn_float_dir, pn_float_name, po_system_configuration,dacname )

%
% Annie Wong, 14 June 2011
% Breck Owens, October 2006
%--------------------------------------------------------------------------

%pn_float_dir='uw/';
%pn_float_name='R5902134';
%po_system_configuration = load_configuration( 'ow_config.txt' );


close all
addpath('/home1/homedir5/perso/ccabanes/dvlpRD/Argo/Lib/ObsInSitu/General/')
addpath('/home1/homedir5/perso/ccabanes/dvlpRD/Argo/Lib/ObsInSitu/General/Coriolis/')
addpath('/home1/homedir5/perso/ccabanes/matlab/lib_downlaod/')

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
LONG = Co.longitude.data;figure
Topo_ficin='/home5/pharos/argo/DMARGO/OW/TOPO/topo.onetenthdeg.nc';
Topo=read_netcdf_allthefile(Topo_ficin);
hold on
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

figure
subplot(2,1,1)


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

minx = max(min(min(x),min(LONG))-5,-100);
maxx = max(max(x),max(LONG))+5;
miny = min(min(y),min(LAT))-5;
maxy=  max(max(y),max(LAT))+5;

if maxy<75; maxy=maxy;end

if  maxy>=75; maxy=75;end


axis([minx maxx miny maxy])

%legend('float','historical points','Location','Best')
%plot(coastdata_x,coastdata_y,'k.-');

for i=1:n
  h=plot(LONG(i),LAT(i),'+','MarkerSize',4);
  set(h,'color',c(i,:));
 % j=text(LONG(i),LAT(i),int2str(PROFILE_NO(i)));
  %set(j,'color',c(i,:),'fontsize',12,'hor','cen');
end

 
set(gca,'FontSize',11)
xlabel('Longitude');
ylabel('Latitude');

filename = ['/home4/begmeil/coriolis_diffusion/ftp/co0508/dac/' strtrim(dacname) '/' strtrim(title_floatname) '/' strtrim(title_floatname) '_prof.nc'];
% disp([num2str(unum(k)) ' ' F.wmo_inst_type.data(1,:)])

F = read_netcdf_allthefile(filename);

title( strcat(title_floatname, ' ARGO +CTD REF : ', po_system_configuration.CONFIGURATION_FILE(1:end-4), ',  (', strtrim(F.pi_name.data(1,:)), '---', strtrim(F.data_centre.data(1,:)),')' ),'interpreter','none' );




% calibration curve (figure 3) --------------------------

if(isempty(find(isnan(cal_SAL)==0))==0) % if no cal exists, terminate here, no plot will be produced

Soffset=cal_SAL-SAL;
avg_Soffset=NaN.*ones(1,n);
avg_Soffset_err=NaN.*ones(1,n);
Staoffset=sta_SAL-SAL;
avg_Staoffset=NaN.*ones(1,n);
avg_Staoffset_err=NaN.*ones(1,n);
avg_mapped_err=NaN.*ones(1,n);

[tlevels, plevels, index, var_s_Thetalevels, Thetalevels] = find_10thetas( SAL, PTMP, PRES, la_ptmp, use_theta_gt, use_theta_lt, use_pres_gt, use_pres_lt, use_percent_gt);

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
   ii=find(isnan(Soffset(:,i))==0);
   if ~isempty(ii)
      
       avg_mapped_err(i)=meanoutnan(mapsalerrors(ii,i));
   else
      
       avg_mapped_err(i) = NaN;
   end
end

% plot salinity time series on theta levels with the smallest S variance (figure 6) ------------

[tlevels, plevels, index, var_s_Thetalevels, Thetalevels] = find_10thetas( SAL, PTMP, PRES, la_ptmp, use_theta_gt, use_theta_lt, use_pres_gt, use_pres_lt, use_percent_gt);

tplot=[1:10];
p=length(tplot);
Sint=NaN.*ones(p,n);
Smap=NaN.*ones(p,n);
Smaperr=NaN.*ones(p,n);
Scal=NaN.*ones(p,n);
Scalerr=NaN.*ones(p,n);
Thetalevel_indexes=NaN.*ones(p,n);

trimPRES=PRES;  % use only manually specified THETA & PRES range ---
trimSAL=SAL;
trimPTMP=PTMP;
trim_mapped_sal=mapped_sal;
trim_mapsalerrors=mapsalerrors;
trim_cal_SAL=cal_SAL;
trim_cal_SAL_err=cal_SAL_err;

jj=find(isnan(la_ptmp)==1);
trimPRES(jj)=NaN;
trimSAL(jj)=NaN;
trimPTMP(jj)=NaN;
trim_mapped_sal(jj)=NaN;
trim_mapsalerrors(jj)=NaN;
trim_cal_SAL(jj)=NaN;
trim_cal_SAL_err(jj)=NaN;

if( isempty(use_theta_lt)==0 & isempty(use_theta_gt)==1 )
  jj=find(trimPTMP>use_theta_lt);
  trimPRES(jj)=NaN;
  trimSAL(jj)=NaN;
  trimPTMP(jj)=NaN;
  trim_mapped_sal(jj)=NaN;
  trim_mapsalerrors(jj)=NaN;
  trim_cal_SAL(jj)=NaN;
  trim_cal_SAL_err(jj)=NaN;
end

if( isempty(use_theta_gt)==0 & isempty(use_theta_lt)==1 )
  jj=find(trimPTMP<use_theta_gt);
  trimPRES(jj)=NaN;
  trimSAL(jj)=NaN;
  trimPTMP(jj)=NaN;
  trim_mapped_sal(jj)=NaN;
  trim_mapsalerrors(jj)=NaN;
  trim_cal_SAL(jj)=NaN;
  trim_cal_SAL_err(jj)=NaN;
end

if( isempty(use_theta_gt)==0 & isempty(use_theta_lt)==0 )
  if(use_theta_gt>use_theta_lt) %the middle band is excluded
    jj=find(trimPTMP<use_theta_gt&trimPTMP>use_theta_lt);
  else
    jj=find(trimPTMP<use_theta_gt|trimPTMP>use_theta_lt);
  end
  trimPRES(jj)=NaN;
  trimSAL(jj)=NaN;
  trimPTMP(jj)=NaN;
  trim_mapped_sal(jj)=NaN;
  trim_mapsalerrors(jj)=NaN;
  trim_cal_SAL(jj)=NaN;
  trim_cal_SAL_err(jj)=NaN;
end

if( isempty(use_pres_lt)==0 & isempty(use_pres_gt)==1 )
  jj=find(trimPRES>use_pres_lt);
  trimPRES(jj)=NaN;
  trimSAL(jj)=NaN;
  trimPTMP(jj)=NaN;
  trim_mapped_sal(jj)=NaN;
  trim_mapsalerrors(jj)=NaN;
  trim_cal_SAL(jj)=NaN;
  trim_cal_SAL_err(jj)=NaN;
end

if( isempty(use_pres_gt)==0 & isempty(use_pres_lt)==1 )
  jj=find(trimPRES<use_pres_gt);
  trimPRES(jj)=NaN;
  trimSAL(jj)=NaN;
  trimPTMP(jj)=NaN;
  trim_mapped_sal(jj)=NaN;
  trim_mapsalerrors(jj)=NaN;
  trim_cal_SAL(jj)=NaN;
  trim_cal_SAL_err(jj)=NaN;
end

if( isempty(use_pres_gt)==0 & isempty(use_pres_lt)==0 )
  if(use_pres_gt>use_pres_lt) %he middle band is excluded
    jj=find(trimPRES<use_pres_gt&trimPRES>use_pres_lt);
  else
    jj=find(trimPRES<use_pres_gt|trimPRES>use_pres_lt);
  end
  trimPRES(jj)=NaN;
  trimSAL(jj)=NaN;
  trimPTMP(jj)=NaN;
  trim_mapped_sal(jj)=NaN;
  trim_mapsalerrors(jj)=NaN;
  trim_cal_SAL(jj)=NaN;
  trim_cal_SAL_err(jj)=NaN;
end

for i=1:n
  for j=tplot
   if(tlevels(j)<max(trimPTMP(:,i))&tlevels(j)>min(trimPTMP(:,i)))
    diffTheta = abs(trimPTMP(:,i)-tlevels(j));
    if isempty(find(~isnan(diffTheta)))
      Thetalevel_indexes(j,i) = NaN;
    else
      Thetalevel_indexes(j,i) = min(find(diffTheta==min(diffTheta)));
    end
   end
  end
end

for i=tplot % build the S matrix for plotting
  for j=1:n
   ti=Thetalevel_indexes(i,j);
   if ~isnan(ti)
     interval=max(ti-1,1):min(ti+1,m); %interval is one above and one below ti
     a = trimPTMP(ti,j) - trimPTMP(interval, j);
     if( trimPTMP(ti,j)>tlevels(i) )
        gg=find(a>0);
        if( ~isempty(gg) )
           b=find(a==min(a(gg))); %find the level with min +ve diff
           ki=interval(b);
        else
           ki=ti;
        end
     end
     if( trimPTMP(ti,j)<tlevels(i) )
        gg=find(a<0);
        if( ~isempty(gg) )
           b=find(-a==min(-a(gg))); %find the level with min -ve diff
           ki=interval(b);
        else
           ki=ti;
        end
     end
     if( trimPTMP(ti,j)==tlevels(i) )
        ki=ti;
     end
     if( ki~=ti&~isnan(trimSAL(ti,j))&~isnan(trimSAL(ki,j))&~isnan(trimPTMP(ti,j))&~isnan(trimPTMP(ki,j)) )
       Sint(i,j) = interp1( [trimPTMP(ti,j), trimPTMP(ki,j)], [trimSAL(ti,j), trimSAL(ki,j)], tlevels(i) );
     else
       Sint(i,j) = trimSAL(ti,j); % interpolate if possible because that is more accurate than using closest points
     end
     if( ki~=ti&~isnan(trim_mapped_sal(ti,j))&~isnan(trim_mapped_sal(ki,j))&~isnan(trimPTMP(ti,j))&~isnan(trimPTMP(ki,j)) )
       Smap(i,j) = interp1( [trimPTMP(ti,j), trimPTMP(ki,j)], [trim_mapped_sal(ti,j), trim_mapped_sal(ki,j)], tlevels(i) );
       Smaperr(i,j) = interp1( [trimPTMP(ti,j), trimPTMP(ki,j)], [trim_mapsalerrors(ti,j), trim_mapsalerrors(ki,j)], tlevels(i) );
     else
       Smap(i,j)=trim_mapped_sal(ti,j); % interpolate if possible because that is more accurate than using closest points
       Smaperr(i,j)=trim_mapsalerrors(ti,j); % interpolate if possible because that is more accurate than using closest points
     end
     if( ki~=ti&~isnan(trim_cal_SAL(ti,j))&~isnan(trim_cal_SAL(ki,j))&~isnan(trimPTMP(ti,j))&~isnan(trimPTMP(ki,j)) )
       Scal(i,j) = interp1( [trimPTMP(ti,j), trimPTMP(ki,j)], [trim_cal_SAL(ti,j), trim_cal_SAL(ki,j)], tlevels(i) );
       Scalerr(i,j) = interp1( [trimPTMP(ti,j), trimPTMP(ki,j)], [trim_cal_SAL_err(ti,j), trim_cal_SAL_err(ki,j)], tlevels(i) );
     else
       Scal(i,j)=trim_cal_SAL(ti,j); % interpolate if possible because that is more accurate than using closest points
       Scalerr(i,j)=trim_cal_SAL_err(ti,j); % interpolate if possible because that is more accurate than using closest points
     end
   end
  end
end
for j=1:n
    ii=find(~isnan(Smaperr(:,j)));
    cov = build_ptmp_cov( tlevels(ii));
    mean_err_mapp(j)=sum(sum(diag(Smaperr(ii,j))*cov))/(length(ii)*length(ii));
end

% Retrouve la correction DM

%keyboard

F = replace_fill_bynan(F);


findDmode = find(F.data_mode.data=='D');
isDmode = F.data_mode.data=='D';

if sum(isDmode)>0
    
    F = extract_profile(F,'N_PROF',findDmode);
    
    F = format_flags_char2num(F);
    
    N_PROF_this_file = size(F.cycle_number.data,1);
    
    
    diff=F.pres_adjusted.data-F.pres.data;
    
    if mean(meanoutnan(diff))~=0
        ij=find(diff~=0&~isnan(diff));
        % compute conductivity from the salinity and raw pressure
        cndr = sw_cndr(F.psal.data(ij),F.temp.data(ij),F.pres.data(ij));
        PRESINI=F.pres_adjusted.data;
        
        
        % recompute the salinity from conductivity and adjusted pressure
        sal = sw_salt(cndr,F.temp.data(ij),F.pres_adjusted.data(ij));
        sal(isnan(sal))=NaN;
        SALINI(ij)=sal;
    end
end
param='psal';
paramad='psal_adjusted';

if sum(isDmode)>0
mean_diff = NaN*zeros(N_PROF_this_file,1);
isnanprof=isnan(F.(param).data)|isnan(F.(paramad).data);
F.(paramad).data(isnanprof)=NaN;
F.(param).data(isnanprof)=NaN;
mean_diff = meanoutnan(F.(paramad).data,2)-meanoutnan(F.(param).data,2);
end
thedate_1cycle=datestr(F.juld.data(1)+datenum('01011950','ddmmyyyy'));
    
    
%keyboard    
subplot(2,1,2)
plot(PROFILE_NO, avg_Soffset, 'b-');
hold on
plot(PROFILE_NO, avg_Soffset, 'g-');
% Plot station by station fit
ok = find(isfinite(avg_Staoffset));
plot(PROFILE_NO(ok), avg_Staoffset(ok), 'r-+','LineWidth',1.5);
if sum(isDmode)>0
plot(F.cycle_number.data,mean_diff,'m-','LineWidth',2)
end
legend('2 x cal error','1 x cal error',['1-1 profile fit (' strrep(po_system_configuration.CONFIGURATION_FILE(1:end-4),'_',' ') ')' ], 'Old DM correction', 'Location', 'SouthOutside','orientation','horizontal');

errorbar(PROFILE_NO, avg_Soffset, 2*avg_Soffset_err,'b')
errorbar(PROFILE_NO, avg_Soffset, avg_Soffset_err,'go-')
%plot(PROFILE_NO(ok), avg_Staoffset(ok), 'r-');

%keyboard
addpath('/home1/corsen/perso/ccabanes/matlab/lib_downlaod/')
coly=[1 0 0];
colerr=[1 0.8 0.8];

myebcc(PROFILE_NO(ok),avg_Staoffset(ok),mean_err_mapp(ok),coly,colerr)

%axis([ 0, max(PROFILE_NO)+1, min([avg_Soffset-avg_Soffset_err,0])-.02, max([avg_Soffset+avg_Soffset_err,0])+.02 ])
%axis([ 0, max(PROFILE_NO)+1, min(avg_Staoffset), max(avg_Staoffset)])
axis tight
plot( [0, max(PROFILE_NO)+1], [0,0], 'k-')
set(gca,'FontSize',10)
xlabel(['float profile number (date of the 1 cycle: ' thedate_1cycle ')' ]);
ylabel('\Delta S')
title( [' vertically-averaged salinity (PSS-78) additive correction \Delta S (red: config', strrep(po_system_configuration.CONFIGURATION_FILE(1:end-4),'_',' '), ')'] );
grid on
if sum(isDmode)>0
plot(F.cycle_number.data,mean_diff,'m-','LineWidth',2)
end

end %if(isempty(find(isnan(cal_SAL)==0))==0) ---------------


drawnow
set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.75,8,9.5]);
print('-depsc ', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_9.eps'));

%  set(gcf,'papertype','usletter','paperunits','inches','paperorientation','portrait','paperposition',[.25,.75,8,9.5]);
%  print('-depsc ', strcat(po_system_configuration.FLOAT_PLOTS_DIRECTORY, pn_float_dir, pn_float_name, '_9.eps'));
%  disp(['!ps2pdf ' po_system_configuration.FLOAT_PLOTS_DIRECTORY pn_float_dir pn_float_name '_9.eps'])
%  eval(['!ps2pdf ' po_system_configuration.FLOAT_PLOTS_DIRECTORY pn_float_dir pn_float_name '_9.eps ' po_system_configuration.FLOAT_PLOTS_DIRECTORY pn_float_dir pn_float_name '_9.pdf'  ])


disp(' ')
disp(' ')



end %if(isempty(find(isnan(mapped_sal)==0))==0) -----------

end %if(isempty(find(isnan(PRES)==0))==0) ------------------


