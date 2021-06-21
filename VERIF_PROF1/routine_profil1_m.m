% -========================================================
%   USAGE :   routine_profil1(floatname,dacname,profnum,PARAM,CONFIG,CAMPAIGN)
%   PURPOSE : plot data from a given float
% -----------------------------------
%   INPUT :
%    floatname  (char)  e.g. '690258'
%    dacname    (char)  e.g.  'coriolis'
%    numcycle   (float array)  e.g. 1 (first cycle 001) see option DIRECTION to consider descending profiles
%                              
%   PARAM   (structure)   input from verif_profil1.m%                            
%   CONFIG   (structure)  general config variables (config.txt)
%   CAMPAIGN  (structure) config variables specific to a campaign  
% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  :  created: V. Thierry - N. David - 08/2007
%               revised: C. Lagadec 2015
%               revised: ccabanes - 2016
%   CALLED SUBROUTINES:
% ========================================================
function  [msg]=routine_profil1(floatname,dacname,profnum,PARAM,CONFIG,CAMPAIGN)


CONFIG.dir_enregistre = [CONFIG.DIR_PLOT '/verif_profil1/' floatname];

if ~exist(CONFIG.dir_enregistre,'dir')
    mkdir(CONFIG.dir_enregistre);
end


% read the float data
FLOAT_SOURCE_NETCDF = CONFIG.(PARAM.DATAREP);
% if strcmp(PARAM.DATATYPE,'adj') == 1
    % FLOAT_SOURCE_NETCDF = CONFIG.DIR_FTP;
% elseif strcmp(PARAM.DATATYPE,'raw') == 1
    % FLOAT_SOURCE_NETCDF = CONFIG.DIR_FTP;
% end

argopath = FLOAT_SOURCE_NETCDF;
repnc = [argopath dacname '/' floatname '/profiles/'];

pfnu=sprintf('%3.3i',profnum);
if strcmp(PARAM.DIRECTION,'D')
pfnu=[pfnu 'D'];
end

fname=['D' floatname '_' pfnu '.nc'];
if exist([repnc fname]) == 2
   file_name=[repnc fname];
else
   fname=['R' floatname '_' pfnu '.nc'];
   if exist([repnc fname]) == 2
      file_name=[repnc fname];
   else
      disp(['Files ' fname  ' does not exist']);
      disp([repnc]);
       msg=['Profile does not exist'];
   return 
   end
end

disp(file_name)

IncludeDescProf=1;
[file_list] = select_float_files_on_ftp(floatname,dacname,FLOAT_SOURCE_NETCDF,'C',IncludeDescProf);
ideF = find(strcmp(file_list,fname));

% select profiles before/after pfnu (PARAM.NB_PROF)
%keyboard
idbef = max(1,ideF-PARAM.NB_PROF);
idaft = min(length(file_list),ideF+PARAM.NB_PROF);

file_list=file_list(idbef:idaft);

[F,Dim,file_list]=create_multi_from_filelist(floatname,dacname,FLOAT_SOURCE_NETCDF,file_list,'Primary',[]);
F = replace_fill_bynan(F);
F = format_flags_char2num(F);

nocycl = find(strcmp(file_list,fname));

% lat_argo = F.latitude.data(nocycl,:);
% lon_argo = F.longitude.data(nocycl,:);
% juld_argo = F.juld.data(nocycl,:);
% psal_argo = F.psal.data(nocycl,:);
% temp_argo = F.temp.data(nocycl,:);
% pres_argo = F.pres.data(nocycl,:);

% psalqc_argo = F.psal_qc.data(nocycl,:)';
% tempqc_argo = F.temp_qc.data(nocycl,:)';
% presqc_argo = F.pres_qc.data(nocycl,:)';


if strcmp(PARAM.DATATYPE,'adj')
F.psal=F.psal_adjusted;
F.psal_qc=F.psal_adjusted_qc;
F.temp=F.temp_adjusted;
F.temp_qc=F.temp_adjusted;
F.pres=F.pres_adjusted;
F.pres_qc=F.pres_adjusted;
end
if strcmp(PARAM.DATAREP,'DIR_DM_FILES')
thepostfix='_adj';
else
thepostfix='';
end



F.psal.data(F.psal_qc.data>3)=NaN;

F.date.data = datevec(datenum(1950,1,1,0,0,0) + F.juld.data);
F.date_str.data = datestr(datenum(1950,1,1,0,0,0) + F.juld.data,1);

F.ptmp0.data = sw_ptmp(F.psal.data,F.temp.data,F.pres.data,0);
F.ptmp0.long_name='Potential Temperature';
F.ptmp0_qc.data=F.psal_qc.data;
F.ptmp1000.data = sw_ptmp(F.psal.data,F.temp.data,F.pres.data,1000);
F.ptmp1000.long_name='Potential Temperature (ref 1000db)';
F.ptmp1000_qc.data=F.psal_qc.data;
%--------------------------------------------------------------------------
% Look for closest shipboard CTD profile
%--------------------------------------------------------------------------
%keyboard
filename = fullfile( CONFIG.DIR_CAMPAIGN, strcat(CAMPAIGN.CAMPAGNE_MAT, CAMPAIGN.POSTFIX ));
CTD_ALL = load( filename );

%lat_campagne = CTD.LAT;
%lon_campagne = CTD.LONG;

dista1 = andoyer(CTD_ALL.LONG,CTD_ALL.LAT, F.longitude.data(nocycl,:),F.latitude.data(nocycl,:));
[distascend,isort]=sort(dista1);
profil_ok=isort(1);
CTD.latitude.data = CTD_ALL.LAT(profil_ok);
CTD.longitude.data = CTD_ALL.LONG(profil_ok);
CTD.juld.data = CTD_ALL.DATES(profil_ok);
CTD.psal.data = CTD_ALL.SAL(profil_ok,:);
CTD.temp.data = CTD_ALL.TEMP(profil_ok,:);
CTD.pres.data = CTD_ALL.PRES(profil_ok,:);

CTD.ptmp0.data  = sw_ptmp(CTD.psal.data,CTD.temp.data,CTD.pres.data,0);
CTD.ptmp1000.data  = sw_ptmp(CTD.psal.data,CTD.temp.data,CTD.pres.data,1000);
CTD.date.data = datevec(double(datenum(1950,1,1,0,0,0) + CTD.juld.data));
CTD.date_str.data = datestr(double(datenum(1950,1,1,0,0,0) + CTD.juld.data),1);

%---------------
% Interpolation 
%---------------

CTDi_on_pres = interpolation_m(CTD,F,'pres',PARAM); % interpolation of CTD variables on float presssure levels F.pres

CTDi_on_ptmp0 = interpolation_m(CTD,F,'ptmp0',PARAM); % interpolation of CTD variables  on  float theta levels F.ptmp0


%---------------
% Computation of PSAL difference on theta levels:
%---------------
difference_psal_theta = -[CTDi_on_ptmp0.psal.data(nocycl,:)-F.psal.data(nocycl,:)];
difference_pres_theta = -[CTDi_on_ptmp0.pres.data(nocycl,:)-F.pres.data(nocycl,:)];
ip = find(abs(difference_pres_theta)<150);


kk = find(F.pres.data(nocycl,ip)>0);
correction_0 = median(difference_psal_theta(ip(kk)));
disp(['Mean Psal difference on float theta levels (pres>0): ' num2str(correction_0)])
kk = find(F.pres.data(nocycl,ip)>500);
correction_500 = median(difference_psal_theta(ip(kk)));
disp(['Mean Psal difference on float theta levels (pres>500): ' num2str(correction_500)])

kk = find(F.pres.data(nocycl,ip)>1000);
correction_1000 = median(difference_psal_theta(ip(kk)));
disp(['Mean Psal difference on float theta levels (pres>1000): ' num2str(correction_1000)])

kk = find(F.pres.data(nocycl,ip)>2000);
correction_2000 = median(difference_psal_theta(ip(kk)));
disp(['Mean Psal difference on float theta levels (pres>2000): ' num2str(correction_2000)])


% difference_temp_brut = [temp_campagne_interpole_brut - temp_argo(1,:)];
% difference_psal = [psal_campagne_interpole_brut - psal_argo(1,:)];

% % "Ref 0 db"

% % [psal_campagne_interpole_ref0,temp_campagne_interpole_ref0] = interpolation(PARAM.LIM_SURF_DEEP,PARAM.STEP_SURF,PARAM.STEP_DEEP,pres_argo,CTD.pres.data,CTD.ptmp0.data,CTD.psal.data);
% difference_temp_ref0 = [temp_campagne_interpole_ref0 - temp_argo_ref0];

% % "Ref 1000 db"

% % [psal_campagne_interpole_ref1000,temp_campagne_interpole_ref1000] = interpolation(PARAM.LIM_SURF_DEEP,PARAM.STEP_SURF,PARAM.STEP_DEEP,pres_argo,CTD.pres.data,CTD.ptmp1000.data,CTD.psal.data);
% difference_temp_ref1000 = [temp_campagne_interpole_ref1000 - temp_argo_ref1000];


% %--------------------------------------------------------------------------
% % Conversion temporelle
% %--------------------------------------------------------------------------

% datetot    = datenum(1950,1,1,0,0,0) + juld_argo;
% date_argo  = datevec(datetot,1);
% date_argo2 = datestr(datetot,1);

% datetot        = datenum(1950,1,1,0,0,0) + double(juld_campagne);
% date_campagne  = datevec(datetot);
% date_campagne2 = datestr(datetot,1);



MAP_VISU.max_lat = max(F.latitude.data(nocycl,:),CTD.latitude.data);
MAP_VISU.max_lon = max(F.longitude.data(nocycl,:),CTD.longitude.data);
MAP_VISU.min_lat = min(F.latitude.data(nocycl,:),CTD.latitude.data);
MAP_VISU.min_lon = min(F.longitude.data(nocycl,:),CTD.longitude.data);
MAP_VISU.zone_visu = [MAP_VISU.min_lat-2 MAP_VISU.max_lat+2 MAP_VISU.min_lon-2 MAP_VISU.max_lon+2];
MAP_VISU.reso='LR';
MAP_VISU.proj='mercator';

%--------------------------------------------------------------------------
%  FIGURE 1 : MAP
%--------------------------------------------------------------------------

hf=figure;
hf.Position=[.25,.75,9.5,8];
hold on
[hf,ha]=fct_pltmap(hf,MAP_VISU.zone_visu,MAP_VISU.reso,MAP_VISU.proj);
set(gca,'fontsize',18);
m_plot(F.longitude.data(nocycl,:),F.latitude.data(nocycl,:),'g.','MarkerSize',20);
m_plot(CTD.longitude.data,CTD.latitude.data,'m.','MarkerSize',20);
title({['Float ' floatname ' cycle ' num2str(nocycl)   ' (' F.date_str.data(nocycl,:) ')'];['vs  CTD from ' CAMPAIGN.camp_name '(' CTD.date_str.data ')']},'FontWeight','bold');
xlabel('Longitude');
ylabel('Latitude');

dista=andoyer(F.longitude.data(nocycl,:),F.latitude.data(nocycl,:),CTD.longitude.data,CTD.latitude.data);
m_text(MAP_VISU.min_lon-1,MAP_VISU.min_lat-1,['Dist: ' num2str(round(dista)) ' km'],'HorizontalAlignment','left','VerticalAlignment','bottom','fontsize',14)
set(gcf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]);
eval(['print -depsc2 ' [CONFIG.dir_enregistre '/' floatname '_map_prof' pfnu  thepostfix '.eps']]);
eval(['print -dpng   ' [CONFIG.dir_enregistre '/' floatname '_map_prof' pfnu   thepostfix '.png']]);

format short g
localisation = [F.date.data(nocycl,:),F.latitude.data(nocycl,:),F.longitude.data(nocycl,:);CTD.date.data,CTD.latitude.data,CTD.longitude.data];
dlmwrite([CONFIG.dir_enregistre '/' floatname '_localisation_' pfnu   '.out'],localisation,' ');

%--------------------------------------------------------------------------
%  FIGURE 2 : MAP & THETA/S
%--------------------------------------------------------------------------
h5 = figure;
h5.Position=[.25,.75,9.5,8];
% flag color
tabcol=[0 1 0 ; 1 1 0 ; 1 0.5 0 ; 1 0 0 ; 0.5 0.5 0 ; 0.5 0.5 0 ; 0.5 0.5 0 ; 0.5 0.5 0; 0.5 0.5 0 ; 0.65 0.65 0.65]; 

PARAM.ptmpref='ptmp0';
[h5,l,b2]=init_figure(h5,F,CTD,floatname,nocycl,pfnu,tabcol,MAP_VISU,PARAM,CAMPAIGN);


% THETA/S diagram
subplot(1,12,[6:12])
hold on
set(gca,'Fontsize',10);

hold on
kk1=(CTD.pres.data>0);
kk2=(F.pres.data(nocycl,:)>0);

for i=1:length(F.latitude.data)
    kk3=F.pres.data(i,:)>0;
    h=plot(F.psal.data(i,kk3),F.(PARAM.ptmpref).data(i,kk3),'k');
end

plot_profile_with_flag(F,'psal','ptmp0',nocycl,'LineWidth',1)


TheXLim=get(gca,'Xlim');
TheYLim=get(gca,'YLim');


plot(CTD.psal.data(kk1),CTD.(PARAM.ptmpref).data(kk1),'.m')
if strcmp(PARAM.DATATYPE,'adj')
ylabel([F.(PARAM.ptmpref).long_name ' Adjusted'])
xlabel([F.psal.long_name ' Adjusted'])   
else
ylabel(F.(PARAM.ptmpref).long_name)
xlabel(F.psal.long_name)
end
title({['Theta / S diagram '];['(Diff. PSAL on theta levels: ' num2str(correction_0) ')']},'FontWeight','bold','Fontsize',11);

grid on;
box on;
axis([TheXLim TheYLim]);
% c=get(gca,'position');
% b=get(l,'position');
% b(2)=c(2);
% set(l,'position',b);
% subplot(1,3,1);
% set(gca,'position',b2);
set(gcf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]); 
eval(['print  -dpsc2 ' [CONFIG.dir_enregistre '/' floatname '_T_S_bathy_' PARAM.ptmpref thepostfix ]  '.eps']);
eval(['print   -dpng ' [CONFIG.dir_enregistre '/' floatname '_T_S_bathy_' PARAM.ptmpref thepostfix]  '.png']);

%--------------------------------------------------------------------------
%  FIGURE 3 : MAP &  PSAL difference
%--------------------------------------------------------------------------

h6 = figure;
h6.Position(3)=h6.Position(3)+h6.Position(3)*30/100;
h6.Position(4)=h6.Position(4)+h6.Position(4)*10/100;
[h6,l,b2]=init_figure(h6,F,CTD,floatname,nocycl,pfnu,tabcol,MAP_VISU,PARAM,CAMPAIGN,'surroundprof','off');

subplot(1,12,[6:12])
hold on
set(gca,'Fontsize',10);

hold on
plot(difference_psal_theta(ip),F.pres.data(nocycl,ip),'+b')
set(gca,'Ydir','reverse')
grid on
box on
if strcmp(PARAM.DATATYPE,'adj')
xlabel('PSAL_ADJ_argo -PSAL_CTD','interpreter','none')
ylabel( [F.pres.long_name ' Adjusted'])
else
xlabel('PSAL_argo -PSAL_CTD','interpreter','none')
ylabel(F.pres.long_name)
end
title({['Diff. PSAL on theta levels: ' num2str(correction_0)]},'FontWeight','bold','Fontsize',11);

set(gcf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]);
eval(['print  -dpsc2 ' [CONFIG.dir_enregistre '/' floatname '_diffPSAL_' PARAM.ptmpref thepostfix]  '.eps']);
eval(['print  -dpng ' [CONFIG.dir_enregistre '/' floatname '_diffPSAL_' PARAM.ptmpref thepostfix]  '.png']);

%--------------------------------------------------------------------------
%  FIGURE 4 : MAP &  PTMP0, PSAL 
%--------------------------------------------------------------------------

h7 = figure;
h7.Position(3)=h7.Position(3)+h7.Position(3)*30/100;
h7.Position(4)=h7.Position(4)+h7.Position(4)*10/100;
[h7,l,b2]=init_figure(h7,F,CTD,floatname,nocycl,pfnu,tabcol,MAP_VISU,PARAM,CAMPAIGN,'surroundprof','on');


% THETA/S diagram
subplot(1,12,[6:8])
hold on
set(gca,'Fontsize',10);

for i=1:length(F.latitude.data)
    kk3=F.pres.data(i,:)>0;
    h=plot(F.(PARAM.ptmpref).data(i,kk3),F.pres.data(i,kk3),'k');
end
kk1=(CTD.pres.data>0);
plot(CTD.(PARAM.ptmpref).data(kk1),CTD.pres.data(kk1),'.m')

plot_profile_with_flag(F,PARAM.ptmpref,'pres',nocycl,'LineWidth',1)
grid on
box on
if strcmp(PARAM.DATATYPE,'adj')

xlabel([F.(PARAM.ptmpref).long_name ' Adj.'])
ylabel([F.pres.long_name ' Adj.'])
else

xlabel(F.(PARAM.ptmpref).long_name)
ylabel(F.pres.long_name)    
end
subplot(1,12,[10:12])
hold on
set(gca,'Fontsize',10);

for i=1:length(F.latitude.data)
    kk3=F.pres.data(i,:)>0;
    h=plot(F.psal.data(i,kk3),F.pres.data(i,kk3),'k');
end
kk1=(CTD.pres.data>0);
plot(CTD.psal.data(kk1),CTD.pres.data(kk1),'.m')
plot_profile_with_flag(F,'psal','pres',nocycl,'LineWidth',1)
grid on
box on

if strcmp(PARAM.DATATYPE,'adj')
xlabel([F.psal.long_name ' Adj.'])
ylabel([F.pres.long_name ' Adj.'])
else
 xlabel(F.psal.long_name)
ylabel(F.pres.long_name)
end
set(gcf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]);
eval(['print  -dpsc2 ' [CONFIG.dir_enregistre '/' floatname '_profiles_' PARAM.ptmpref thepostfix]  '.eps']);
eval(['print  -dpng ' [CONFIG.dir_enregistre '/' floatname '_profiles_' PARAM.ptmpref thepostfix]  '.png']);

%--------------------------------------------------------------------------
% THETA/S  with ZOOM
%--------------------------------------------------------------------------

h8 = figure;
h8.Position(3)=h8.Position(3)+h8.Position(3)*30/100;
h8.Position(4)=h8.Position(4)+h8.Position(4)*10/100;
%[h8,l,b2]=init_figure(h8,F,CTD,floatname,nocycl,pfnu,tabcol,MAP_VISU,PARAM,CAMPAIGN,'surroundprof','on');

subplot(1,3,1)
set(gca,'Fontsize',10);
hold on
plot(CTD.psal.data,CTD.(PARAM.ptmpref).data,'.m')
nstop=0;
for vqc=1:10       
    iqc=find(F.psal_qc.data(nocycl,:)==vqc);
    if isempty(iqc)==0&nstop==0
    hqc=plot(F.psal.data(nocycl,iqc),F.(PARAM.ptmpref).data(nocycl,iqc),'o');
    set(hqc,'color',tabcol(vqc,:),'markerfacecolor',tabcol(vqc,:),'markersize',5);
    hold on
    nstop=1;
    end
end

TheXLim=get(gca,'Xlim');
TheYLim=get(gca,'YLim');
nstop=0;
ileg=0;
for i=1:length(F.latitude.data)
    kk3=F.pres.data(i,:)>0;
    h=plot(F.psal.data(i,kk3),F.(PARAM.ptmpref).data(i,kk3),'k');
    if sum(kk3)~=0&nstop==0
        ileg=1;
        nstop=1;
    end
end

for vqc=1:10       
    iqc=find(F.psal_qc.data(nocycl,:)==vqc);
	if isempty(iqc)==0
    hqc=plot(F.psal.data(nocycl,iqc),F.(PARAM.ptmpref).data(nocycl,iqc),'o');
    set(hqc,'color',tabcol(vqc,:),'markerfacecolor',tabcol(vqc,:),'markersize',5);
    hold on
	end
end

plot(CTD.psal.data,CTD.(PARAM.ptmpref).data,'.m')
if strcmp(PARAM.DATATYPE,'adj')

ylabel([F.(PARAM.ptmpref).long_name ' Adjusted'])
xlabel([F.psal.long_name ' Adjusted'])
else
 
ylabel(F.(PARAM.ptmpref).long_name)
xlabel(F.psal.long_name)   
end

grid on 
box on
axis([TheXLim TheYLim]);

% Moved here by T. Reynaud 18.09.2020
if ileg==1
    legend([CAMPAIGN.camp_name],['Argo ' pfnu],['Argo profiles +/- ' num2str(PARAM.NB_PROF) ' cycle(s)'] ,'location','SouthOutside'); 
end

title({['Float ' floatname ' cycle ' pfnu   ' (' F.date_str.data(nocycl,:) ')'];['vs  CTD from ' CAMPAIGN.camp_name '(' CTD.date_str.data ')']},'fontsize',11,'FontWeight','bold');

subplot(1,12,[6:12])
hold on
set(gca,'Fontsize',10);

hold on
kk1=(CTD.pres.data>PARAM.ZOOM);
kk2=(F.pres.data(nocycl,:)>PARAM.ZOOM);
%keyboard

plot(CTD.psal.data(kk1),CTD.(PARAM.ptmpref).data(kk1),'.m')
%plot(psal_campagne_interpole_brut(kk2),temp_campagne_interpole_compil(ind1,kk2),'k+')


for vqc=1:10  
%    iqc=find(max(flag) == vqc)
    iqc=find(F.psal_qc.data(nocycl,:)==vqc&kk2);
    hqc=plot(F.psal.data(nocycl,iqc),F.(PARAM.ptmpref).data(nocycl,iqc),'o');
    set(hqc,'color',tabcol(vqc,:),'markerfacecolor',tabcol(vqc,:),'markersize',5);
    hold on
end

%legend([camp_name],[camp_name ' interp.'],['Argo ' nocycl] ,'location','NorthEast');
nstop=0;
for i=1:length(F.latitude.data)
    kk3=F.pres.data(i,:)>PARAM.ZOOM;
    h=plot(F.psal.data(i,kk3),F.(PARAM.ptmpref).data(i,kk3),'k');

end
%plot(psal_campagne_interpole_brut(kk2),temp_campagne_interpole_compil(ind1,kk2),'k+')

for vqc=1:10  
%    iqc=find(max(flag) == vqc)
    iqc=find(F.psal_qc.data(nocycl,:)==vqc&kk2);
    hqc=plot(F.psal.data(nocycl,iqc),F.(PARAM.ptmpref).data(nocycl,iqc),'o');
    set(hqc,'color',tabcol(vqc,:),'markerfacecolor',tabcol(vqc,:),'markersize',5);
    hold on
end

plot(CTD.psal.data(kk1),CTD.(PARAM.ptmpref).data(kk1),'.m');

% Moved here by T. Reynaud 18.06.2020
TheXLim=get(gca,'Xlim');
TheYLim=get(gca,'YLim');

%ylim([1 6])
if strcmp(PARAM.DATATYPE,'adj')
xlabel([F.(PARAM.ptmpref).long_name ' Adjusted'])
ylabel([F.psal.long_name ' Adjusted'])
else
ylabel(F.(PARAM.ptmpref).long_name)
xlabel(F.psal.long_name)
end
grid on
box on
axis([TheXLim TheYLim])
title({['Zoom on depths deeper than ' num2str(PARAM.ZOOM) 'm'];['Diff. PSAL on theta levels: ' num2str(correction_0)]},'FontWeight','bold','Fontsize',11);

set(gcf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]);
eval(['print  -dpsc2 ' [CONFIG.dir_enregistre '/' floatname '_T_S_zoom_' PARAM.ptmpref thepostfix ]  '.eps']);
eval(['print   -dpng ' [CONFIG.dir_enregistre '/' floatname '_T_S_zoom_' PARAM.ptmpref thepostfix]  '.png']);
    

	
%%%% init_figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h5,legende,b2]=init_figure(h5,F,CTD,floatname,nocycl,pfnu,tabcol,MAP_VISU,PARAM,CAMPAIGN,varargin)
n=length(varargin);

if n/2~=floor(n/2)
    error('check the imput arguments')
end

f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);
surroundprof='on';
if isfield(s,'surroundprof')==1;surroundprof=s.surroundprof;end;

% build the legend:
subplot(1,3,1)
set(gca,'Fontsize',10);
hold on
hjb=plot(CTD.psal.data,CTD.(PARAM.ptmpref).data,'o');
set(hjb,'color','m','markerfacecolor','m','markersize',5);
nstop=0;

for vqc=1:10       
    iqc=find(F.psal_qc.data(nocycl,:)==vqc);
    if isempty(iqc)==0& nstop==0
		hqc = plot(F.psal.data(nocycl,iqc),F.(PARAM.ptmpref).data(nocycl,iqc),'b-+');
		%set(hqc,'color',tabcol(vqc,:),'markerfacecolor',tabcol(vqc,:),'markersize',5);
		%set(hqc,'markerfacecolor',tabcol(vqc,:),'markersize',5);
		hold on
		nstop=1;
		vqc_save=vqc;
    end
end

nstop=0;
for i=1:length(F.latitude.data)
    kk3=F.pres.data(i,:)>0;
    h=plot(F.psal.data(i,kk3),F.(PARAM.ptmpref).data(i,kk3),'k');
    if sum(kk3)~=0&nstop==0        
        if strcmp(surroundprof,'on')
            ileg=0;
            % T. Reynaud 18.09.2020
            %l=legend([CAMPAIGN.camp_name],['Argo ' pfnu],['Argo profiles +/- ' num2str(PARAM.NB_PROF) ' cycle(s)'] ,'location','SouthOutside');
        else
            ileg=1;
            % T. Reynaud 18.09.2020
            %l=legend([CAMPAIGN.camp_name],['Argo ' pfnu] ,'location','SouthOutside');          
        end
        nstop=1;
    end
end   


for vqc=1:10       
    iqc=find(F.psal_qc.data(nocycl,:)==vqc);
	hqc = plot(F.psal.data(nocycl,iqc),F.(PARAM.ptmpref).data(nocycl,iqc),'o');
	set(hqc,'color',tabcol(vqc,:),'markerfacecolor',tabcol(vqc,:),'markersize',5);
	hold on
end

% location MAP:
MAP_VISU.zone_visu = [MAP_VISU.min_lat-1 MAP_VISU.max_lat+1 MAP_VISU.min_lon-1 MAP_VISU.max_lon+1];
[h5,ha]=fct_pltmap(h5,MAP_VISU.zone_visu,MAP_VISU.reso,MAP_VISU.proj);
b2=get(gca,'position');

m_plot(F.longitude.data(nocycl,:),F.latitude.data(nocycl,:),'b.','MarkerSize',20);
m_plot(CTD.longitude.data,CTD.latitude.data,'m.','MarkerSize',20);
title({['Float ' floatname ' cycle ' pfnu   ' (' F.date_str.data(nocycl,:) ')'];['vs  CTD from ' CAMPAIGN.camp_name '(' CTD.date_str.data ')']},'fontsize',11,'FontWeight','bold');
xlabel('Longitude');
ylabel('Latitude');

dista=andoyer(F.longitude.data(nocycl,:),F.latitude.data(nocycl,:),CTD.longitude.data,CTD.latitude.data);
m_text(MAP_VISU.min_lon-0.25,MAP_VISU.min_lat-0.25,['Dist: ' num2str(round(dista)) ' km'],'HorizontalAlignment','left','VerticalAlignment','bottom','fontsize',10,'fontweight','bold')

% Legend moved here by T. Reynaud 18.09.2020
if ileg==0
    %l=legend([CAMPAIGN.camp_name],['Argo ' pfnu],['Argo profiles +/- ' num2str(PARAM.NB_PROF) ' cycle(s)'] ,'location','SouthOutside');
    l=legend(['Argo ' pfnu],[CAMPAIGN.camp_name],['Argo profiles +/- ' num2str(PARAM.NB_PROF) ' cycle(s)'] ,'location','SouthOutside');
elseif ileg==1
    %l=legend([CAMPAIGN.camp_name],['Argo ' pfnu] ,'location','SouthOutside');
    l=legend(['Argo ' pfnu] ,[CAMPAIGN.camp_name],'location','SouthOutside');
end
b=get(l,'position');
legende=l;

