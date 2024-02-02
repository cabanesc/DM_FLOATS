% -========================================================
%   USAGE : PLOT_FLOAT_and_FLOAT(floatname,dacname,floatnameref,dacnameref,campaign_name varargin)
%   PURPOSE : compare two Argo floats (map on theta levels, time series on theta levels theta/S diagram)
% -----------------------------------
%   INPUT :
%    floatname  (char)  e.g. '6902881'
%    dacname    (char) e.g.  'coriolis'
%    floatnameref  (char)  e.g. '6902882'
%    dacnameref    (char) e.g.  'coriolis'
%    campaign_name  (char) e.g 'config_ovide18.txt'
%   OPTIONNAL INPUT :
%     'CYCLE_MIN' (float) minimum float's cycle to consider (default=1)
%     'CYCLE_MAX' (float) maximum float's cycle to consider  (default=last cycle available)
%     'REFCYCLE_MIN' (float) minimum floatref's cycle to consider (default=1)
%     'REFCYCLE_MAX' (float) maximum floareft's cycle to consider  (default=last cycle available)
%     'DATATYPE'  (char)  'raw' if raw PARAM are plotted; 'adj' if PARAM_ADJ are plotted  (default is 'raw
%     'DATAREP'   (char)   'DIR_DM_FILES'
%     'MAX_DIST'  (float)  Maximum distance between two profiles to be compared, default 50km
%     'MIN_DEPTH' (float)  Minimum depth to consider to compute the mean salinity difference, default 1000
%     'PRINT'    (logical) PRINT=1 is figure is saved  PRINT=0 otherwise
%     (default=0)
% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  : created (2020) ccabanes from PLOT_REF_and_FLOAT.m
%
%   CALLED SUBROUTINES:
% ========================================================

function PLOT_FLOAT_and_FLOAT(floatname,dacname,floatnameref,dacnameref,campaign_name, varargin)

close all
%init_path

CONFIG=load_configuration('config.txt');

CONFIG.pathwmobox=[CONFIG.DIR_OWC CONFIG.VERSION_OWC '/data/'];

if exist([CONFIG.pathwmobox 'climatology/historical_ctd/bad_data_point.mat'])
    load([CONFIG.pathwmobox 'climatology/historical_ctd/bad_data_point.mat']);
    CONFIG.ctd_to_exclude=bad_ctd.source;
else
    CONFIG.ctd_to_exclude={''};
end

% To save plots
CONFIG.plotpathini=[CONFIG.DIR_PLOT '/ref_database/'];
CONFIG.plotpath=[CONFIG.plotpathini '/' floatname '/'];
if exist(CONFIG.plotpath)~=7
    mkdir(CONFIG.plotpath)
end

n=length(varargin);

if n/2~=floor(n/2)
    error('check the imput arguments')
end

f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);

NcVar.latitude.name='LATITUDE';
NcVar.longitude.name='LONGITUDE';
NcVar.cycle_number.name='CYCLE_NUMBER';
CONFIG.DATATYPE='raw';
CONFIG.DATAREP='DIR_FTP';
if isfield(s,'DATATYPE')==1;CONFIG.DATATYPE=s.DATATYPE;end;

if strcmp(CONFIG.DATATYPE,'adj') == 1
    isbest_value=2;
elseif strcmp(CONFIG.DATATYPE,'raw') == 1
    isbest_value=0;
end
isbest_value
if isfield(s,'DATAREP')==1;CONFIG.DATAREP=s.DATAREP;end;
if strcmp(CONFIG.DATAREP,'DIR_FTP')==1
    DATAREP=CONFIG.DIR_FTP;
elseif strcmp(CONFIG.DATAREP,'DIR_DM_FILES')==1
    DATAREP=CONFIG.DIR_DM_FILES;
end


%keyboard
IncludeDescProf=1;
[file_list] = select_float_files_on_ftp(floatname,dacname,DATAREP,'C',IncludeDescProf);
[file_list_ref] = select_float_files_on_ftp(floatnameref,dacnameref,DATAREP,'C',IncludeDescProf);

% % select profiles before/after pfnu (PARAM.NB_PROF)
% %keyboard
% idbef = max(1,ideF-PARAM.NB_PROF);
% idaft = min(length(file_list),ideF+PARAM.NB_PROF);

% file_list=file_list(idbef:idaft);

% [F,Dim,file_list]=create_multi_from_filelist(floatname,dacname,FLOAT_SOURCE_NETCDF,file_list,'Primary',[]);
% F = replace_fill_bynan(F);
% F = format_flags_char2num(F);


% F=read_netcdf_allthefile([CONFIG.DIR_FTP dacname '/' floatname '/' floatname '_prof.nc'],NcVar);
% Fref=read_netcdf_allthefile([CONFIG.DIR_FTP dacnameref '/' floatnameref '/' floatnameref '_prof.nc'],NcVar);

% default CONFIG
CONFIG.CYCLE_MIN=1;
CONFIG.CYCLE_MAX=length(file_list);
CONFIG.REFCYCLE_MIN=1;
CONFIG.REFCYCLE_MAX=length(file_list_ref);
CONFIG.PRINT=0;
CONFIG.ZOOM_TPOT=4;
CONFIG.MAX_DIST=50;
CONFIG.MIN_DEPTH=1000;

% Input CONFIG
if isfield(s,'CYCLE_MIN')==1;CONFIG.CYCLE_MIN=s.CYCLE_MIN;end
if isfield(s,'CYCLE_MAX')==1;CONFIG.CYCLE_MAX=s.CYCLE_MAX;end
if isfield(s,'REFCYCLE_MIN')==1;CONFIG.REFCYCLE_MIN=s.REFCYCLE_MIN;end
if isfield(s,'REFCYCLE_MAX')==1;CONFIG.REFCYCLE_MAX=s.REFCYCLE_MAX;end
if isfield(s,'ZOOM_TPOT')==1;CONFIG.ZOOM_TPOT=s.ZOOM_TPOT;end
if isfield(s,'MAX_DIST')==1;CONFIG.MAX_DIST=s.MAX_DIST;end
if isfield(s,'MIN_DEPTH')==1;CONFIG.MIN_DEPTH=s.MIN_DEPTH;end

if isfield(s,'PRINT')==1;CONFIG.PRINT=s.PRINT;end
CONFIG

%file_list=file_list(CONFIG.CYCLE_MIN:CONFIG.CYCLE_MAX);
%file_list_ref=file_list_ref(CONFIG.REFCYCLE_MIN:CONFIG.REFCYCLE_MAX);

[F,Dim,file_list]=create_multi_from_filelist(floatname,dacname,DATAREP,file_list,'Primary',[]);
F = replace_fill_bynan(F);
F = format_flags_char2num(F);

if isbest_value==1
    F =construct_best_param(F ,{'temp','pres','psal'},F);
    F.psal=F.psal_best;F.temp=F.temp_best;F.pres=F.pres_best;
    F.psal.data(F.psal_qc.data>1)=NaN;
end
if isbest_value==2
    %F =construct_best_param(F ,{'temp','pres','psal'},F);
    F.psal=F.psal_adjusted;F.temp=F.temp_adjusted;F.pres=F.pres_adjusted;
    F.psal.data(F.psal_adjusted_qc.data>3)=NaN;
end

[Fref,Dimref,file_list_ref]=create_multi_from_filelist(floatnameref,dacnameref,DATAREP,file_list_ref,'Primary',[]);
Fref = replace_fill_bynan(Fref);
Fref = format_flags_char2num(Fref);


if isbest_value==1
    Fref =construct_best_param(Fref ,{'temp','pres','psal'},F);
    Fref.psal=Fref.psal_best;Fref.temp=Fref.temp_best;Fref.pres=Fref.pres_best;
    Fref.psal.data(Fref.psal_qc.data>1)=NaN;
end
if isbest_value==2
    %F =construct_best_param(F ,{'temp','pres','psal'},F);
    Fref.psal=Fref.psal_adjusted;Fref.temp=Fref.temp_adjusted;Fref.pres=Fref.pres_adjusted;
    Fref.psal.data(Fref.psal_adjusted_qc.data>3)=NaN;
end


% selection des cycles

ideb=find(F.cycle_number.data==CONFIG.CYCLE_MIN);
if length(ideb)>1;    ideb=ideb(1); end   
if isempty(ideb);     ideb=1      ;end
ifin=find(F.cycle_number.data==CONFIG.CYCLE_MAX);
if length(ifin)>1;    ifin=ifin(end)        ; end
if isempty(ifin);     ifin=length(file_list); end
refideb=find(Fref.cycle_number.data==CONFIG.REFCYCLE_MIN);
if length(refideb)>1;    refideb=refideb(1); end
if isempty(refideb);     refideb=1      ;end
refifin=find(Fref.cycle_number.data==CONFIG.REFCYCLE_MAX);
if length(refifin)>1;    refifin=refifin(end)        ; end
if isempty(refifin);     refifin=length(file_list_ref); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DES DONNEES des deux flotteurs

%--------------------------------------------------------------------------
%  FIGURE 1 : MAP et theta/S
%--------------------------------------------------------------------------

% plot diagramme theta/S
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
h1=plot_float_diag(floatnameref,dacnameref,CONFIG,campaign_name,refideb,refifin,isbest_value,'.r');
h2=plot_float_diag(floatname,dacname,CONFIG,campaign_name,ideb,ifin,isbest_value,'.b');

legend([h1;h2],{[floatnameref ' cycles :' num2str(Fref.cycle_number.data(refideb)) '-' num2str(Fref.cycle_number.data(refifin))];[floatname ' cycles :' num2str(F.cycle_number.data(ideb)) '-' num2str(F.cycle_number.data(ifin))] },'Location','best')
%keyboard

titplt=[floatname '_and_' floatnameref '_thetaS.pdf'];


subplot(1,2,1)
MAP_VISU.max_lat = max(max(F.latitude.data(ideb:ifin)),max(Fref.latitude.data(refideb:refifin)));
MAP_VISU.max_lon = max(max(F.longitude.data(ideb:ifin)),max(Fref.longitude.data(refideb:refifin)));
MAP_VISU.min_lat = min(min(F.latitude.data(ideb:ifin)),min(Fref.latitude.data(refideb:refifin)));
MAP_VISU.min_lon = min(min(F.longitude.data(ideb:ifin)),min(Fref.longitude.data(refideb:refifin)));
MAP_VISU.zone_visu = [MAP_VISU.min_lat-2 MAP_VISU.max_lat+2 MAP_VISU.min_lon-2 MAP_VISU.max_lon+2];
MAP_VISU.reso='LR';
MAP_VISU.proj='mercator';

% MAP
%%%%%
hf=gca;
hold on
fct_pltmap_new(MAP_VISU.zone_visu,MAP_VISU.reso,MAP_VISU.proj);
%set(gca,'fontsize',18);
m_plot(F.longitude.data(ideb:ifin),F.latitude.data(ideb:ifin),'.b','MarkerSize',20);
m_plot(Fref.longitude.data(refideb:refifin),Fref.latitude.data(refideb:refifin),'.r','MarkerSize',10);
%xlabel('Longitude');
%ylabel('Latitude');

if CONFIG.PRINT==1
    eval(['print -dpdf ' CONFIG.plotpath  titplt]);
end

% recherche des profils les plus proches.
disp('recherche des profils les plus proches.')
iprofref=NaN*zeros(1,length(F.longitude.data));
vminv=NaN*zeros(1,length(F.longitude.data));

% calcule les distances entre chaque profil des deux flotteurs
% et pour chaque profil du flotteur analyse donne le profil de reference le
% plus proche => iprofref et la distance correspondante => vminv
%for iprof=1:length(F.longitude.data)
for iprof=ideb:ifin
    %figure
    tabdist=andoyer(Fref.longitude.data(refideb:refifin),Fref.latitude.data(refideb:refifin),F.longitude.data(iprof),F.latitude.data(iprof));
    [vmin,io]=sort(tabdist);
    iprofref(iprof)=io(1); %  numero du profil de reference le plus proche du profil iprof
    vminv(iprof)=vmin(1);  %  distance entre le profil iprof et  le plus proche du flotteur de reference
end

clear profan_tokeep distanref
% un meme profil du flotteur de reference peut etre selectionné pour
% plusieurs profils du flotteur analyse => profref_tokeep
% Pour chaque profil de reference selectionne, recherche les profils du
% flotteur analyse les plus proche => profan_tokeep et la distance avec le
% profil de reference =>distanref
%

[profref_tokeep,iref_tokeep]=unique(iprofref);
profref_tokeep(isnan(profref_tokeep))=[];
iref_tokeep(isnan(profref_tokeep))=[];
for ik=1:length(profref_tokeep)                % pour chacun on ne garde qu'un profil du flotteur analysé : le plus proche
    ii=find(iprofref==profref_tokeep(ik));     % ii index des profils du flotteur analyse
    [dmin,dd]=sort(vminv(ii));
    profan_tokeep(ik)=ii(dd(1));
    distanref(ik)=vminv(ii(dd(1)));             % doit etre egal a vminv(profan_tokeep(ik)) distance entre profil analyse et profil de ref
end

% critere de distance maximale entre profil proche
profan_tokeep(distanref>CONFIG.MAX_DIST)=[];
profref_tokeep(distanref>CONFIG.MAX_DIST)=[];
distanref(distanref>CONFIG.MAX_DIST)=[];


% elimine les cycles qui ne sont pas selectionnes au depart
to_elim = profan_tokeep>ifin|profan_tokeep<ideb|profref_tokeep>refifin|profref_tokeep<refideb;
profan_tokeep(to_elim)=[];
profref_tokeep(to_elim)=[];
distanref(to_elim)=[];

MAP_VISU.max_lat = max(max(F.latitude.data),max(Fref.latitude.data));
MAP_VISU.max_lon = max(max(F.longitude.data),max(Fref.longitude.data));
MAP_VISU.min_lat = min(min(F.latitude.data),min(Fref.latitude.data));
MAP_VISU.min_lon = min(min(F.longitude.data),min(Fref.longitude.data));
MAP_VISU.zone_visu = [MAP_VISU.min_lat-2 MAP_VISU.max_lat+2 MAP_VISU.min_lon-2 MAP_VISU.max_lon+2];
MAP_VISU.reso='LR';
MAP_VISU.proj='mercator';

scale=(MAP_VISU.max_lon-MAP_VISU.min_lon)/100;

% MAP
%%%%%
figure(2)
hf=gca;
hold on
fct_pltmap_new(MAP_VISU.zone_visu,MAP_VISU.reso,MAP_VISU.proj);
%set(gca,'fontsize',18);
m_plot(F.longitude.data,F.latitude.data,':b');
m_plot(Fref.longitude.data,Fref.latitude.data,':r');
% j=m_text(F.longitude.data(1)+scale,F.latitude.data(1),int2str(F.cycle_number.data(1)));
% set(j,'color','b','fontsize',10,'hor','cen');
% j=m_text(F.longitude.data(end)+scale,F.latitude.data(end),int2str(F.cycle_number.data(end)));
% set(j,'color','b','fontsize',10,'hor','cen');
% j=m_text(Fref.longitude.data(1)+scale,Fref.latitude.data(1),int2str(Fref.cycle_number.data(1)));
% set(j,'color','r','fontsize',10,'hor','cen');
% j=m_text(Fref.longitude.data(end)+scale,Fref.latitude.data(end),int2str(Fref.cycle_number.data(end)));
% set(j,'color','r','fontsize',10,'hor','cen');
for k=1:10:length(F.longitude.data)
   j=m_text(F.longitude.data(k)+scale,F.latitude.data(k),int2str(F.cycle_number.data(k)));
   set(j,'color','b','fontsize',8,'hor','cen'); 
end
for k=1:10:length(Fref.longitude.data)
   j=m_text(Fref.longitude.data(k)+scale,Fref.latitude.data(k),int2str(Fref.cycle_number.data(k)));
   set(j,'color','r','fontsize',8,'hor','cen'); 
end
m_scatter(F.longitude.data(ideb:ifin),F.latitude.data(ideb:ifin),5,'ob');
m_scatter(Fref.longitude.data(refideb:refifin),Fref.latitude.data(refideb:refifin),5,'or');

m_scatter(F.longitude.data(profan_tokeep),F.latitude.data(profan_tokeep),30,'ob','filled');
m_scatter(Fref.longitude.data(profref_tokeep),Fref.latitude.data(profref_tokeep),30,'or','filled');
title({['Compared profiles (Blue - ' floatname ' and Red - ' floatnameref ')' ];['Maximum distance between profiles :' num2str(CONFIG.MAX_DIST) 'km']})

%--------------------------------------------------------------------------
% FIGURE 3 sera tracé après
figure(3)


%--------------------------------------------------------------------------
%  FIGURE 4 et suivantes : Difference PSAL on theta levels
%--------------------------------------------------------------------------

pres_med=[1:100:4000];
diff_med=zeros(length(pres_med)-1,length([ideb:ifin]));
themean=NaN*zeros(1,length(profan_tokeep));



for ifloat=1:length(profan_tokeep)
    figure
    t1 = tiledlayout(1,3);
    iprof=profan_tokeep(ifloat);
    iov=profref_tokeep(ifloat);
    % calcul de la difference psal sur niveaux theta
    tpot=sw_ptmp(F.psal.data(iprof,:),F.temp.data(iprof,:),F.pres.data(iprof,:),0);
    sig0=sw_pden(F.psal.data(iprof,:),F.temp.data(iprof,:),F.pres.data(iprof,:),0)-1000;
    psal=F.psal.data(iprof,:);
    Fref.tpot.data(iov,:)=sw_ptmp(Fref.psal.data(iov,:),Fref.temp.data(iov,:),Fref.pres.data(iov,:),0);
    [S_h,P_h]=interp_climatology(Fref.psal.data(iov,:)',Fref.tpot.data(iov,:)',Fref.pres.data(iov,:)',psal,tpot,F.pres.data(iprof,:));
    
    difference_psal_theta=[F.psal.data(iprof,:)-S_h'];
    
    difference_pres_theta=[F.pres.data(iprof,:)-P_h'];
    for ik=1:length(pres_med)-1
        diff_med(ik,iprof)=mean(difference_psal_theta(F.pres.data(iprof,:)>pres_med(ik)&F.pres.data(iprof,:)<pres_med(ik+1)));
    end
    %keyboard
    ip=find(abs(difference_pres_theta)<1000);
    kk=find(F.pres.data(iprof,ip)>CONFIG.MIN_DEPTH );
    if length(ip(kk))>1
        themean(ifloat)=mean(difference_psal_theta(ip(kk)));
    else
        themean(ifloat)=NaN;
    end
    nexttile(t1)
    %subplot(1,2,1)
    hold on
     plot(F.temp.data(iprof,:),F.pres.data(iprof,:),'b','LineWidth',2)
     plot(Fref.temp.data(iov,:),Fref.pres.data(iov,:),'r','LineWidth',2)
    set(gca,'Ydir','reverse')
    
    grid on
    box on
    xlabel(F.temp.name)
    ylabel(F.pres.name)
    
    nexttile(t1)
    %subplot(1,2,1)
    hold on
%     plot(F.psal.data(iprof,:),F.pres.data(iprof,:),'b','LineWidth',2)
%     plot(Fref.psal.data(iov,:),Fref.pres.data(iov,:),'r','LineWidth',2)
    plot(F.psal.data(iprof,:),tpot,'b','LineWidth',2)
    plot(Fref.psal.data(iov,:),Fref.tpot.data(iov,:),'r','LineWidth',2)
    %set(gca,'Ydir','reverse')
    grid on
    box on
    xlabel(F.psal.name)
    %ylabel(F.pres.long_name)
    ylabel('Potential Temperature')
    nexttile(t1)
    %subplot(1,2,2)
    hold on
    plot(difference_psal_theta(ip),F.pres.data(iprof,ip),'k','LineWidth',2)
    set(gca,'Ydir','reverse')
    grid on
    box on
    xlabel('Diff PSAL on theta levels')
    ylabel(F.pres.long_name)
    % subplot(1,3,3)
    % plot(difference_pres_theta,F.pres.data(iprof,:),'k','LineWidth',2)
    % set(gca,'Ydir','reverse')
    % grid on
    % box on
    % xlabel('Diff PRES on theta levels')
    % ylabel(F.pres.long_name)
    title(t1,{['Blue - ' floatname ' (cycle ' num2str(F.cycle_number.data(iprof)) strtrim(F.direction.data(iprof)) ') and - Red - ' floatnameref ' (cycle ' num2str(Fref.cycle_number.data(iov)) strtrim(Fref.direction.data(iov)) ').'];['  Mean diff <' num2str(CONFIG.MIN_DEPTH) ' db : ' num2str(themean(ifloat))];['Distance between profiles : ' num2str(distanref(ifloat)) 'km'] })
    
end


figure(3)
hold on
p=plot(F.cycle_number.data(profan_tokeep),themean,'+')
xlim = get(gca,'Xaxis')
xlim.Limits=[F.cycle_number.data(1) F.cycle_number.data(end)+5];
p1=plot(F.cycle_number.data(profan_tokeep),repmat(mean(themean),1,length(profan_tokeep)),'--g','linewidth',2)
%plot(meanoutnan(diff_med'),(pres_med(1:end-1)+pres_med(2:end))/2)
%set(gca,'Ydir','reverse')
grid on
box on
if length(profan_tokeep)>0
legend(p1,['mean : ' num2str(mean(themean))],'Location','SouthOutside')
else
    text(1,0.5,['   No close profiles within ' num2str(CONFIG.MAX_DIST) 'km'])
end
    
title ({['Salinity difference below ' num2str(CONFIG.MIN_DEPTH)]; ['between the closest profiles from : ' floatname ' - ' floatnameref]})
ylabel('\Delta PSAL')
xlabel([floatname ': Cycle Number'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hl=plot_float_diag(floatname,dacname,CONFIG,campaign_name,ideb,ifin,Isbest,Marker)

vertical_sampling_scheme='Primary sampling';
IncludeDescProf=1;

if strcmp(CONFIG.DATAREP,'DIR_FTP')==1
    DATAREP=CONFIG.DIR_FTP;
elseif strcmp(CONFIG.DATAREP,'DIR_DM_FILES')==1
    DATAREP=CONFIG.DIR_DM_FILES;
end

[file_list] = select_float_files_on_ftp(floatname,dacname,DATAREP,'C',IncludeDescProf);
[F,Dim,thelist_ext2] = create_multi_from_filelist(floatname,dacname,DATAREP,file_list,vertical_sampling_scheme,'');

% filename=[CONFIG.DIR_FTP  dacname '/' floatname '/' floatname '_prof.nc'];
% F=read_netcdf_allthefile(filename);
F = replace_fill_bynan(F);
F = format_flags_char2num(F);
F.psal.data(F.psal_qc.data>3)=NaN;
if Isbest==1
    F =construct_best_param(F ,{'temp','pres','psal'},F);
    F.psal=F.psal_best;F.temp=F.temp_best;F.pres=F.pres_best;
    F.psal.data(F.psal_qc.data>1)=NaN;
end
if Isbest==2
    %F =construct_best_param(F ,{'temp','pres','psal'},F);
    F.psal=F.psal_adjusted;F.temp=F.temp_adjusted;F.pres=F.pres_adjusted;
    F.psal.data(F.psal_adjusted_qc.data>3)=NaN;
end
F.tpot=F.temp;
F.tpot.data = sw_ptmp(F.psal.data,F.temp.data,F.pres.data,0);

subplot(1,2,2)
plot(F.psal.data(ideb:ifin,:)',F.tpot.data(ideb:ifin,:)',Marker);
hold on
hl=plot(F.psal.data(ideb,:)',F.tpot.data(ideb,:)',Marker);
ylim=get(gca,'Ylim');

set(gca,'Ylim',[ylim(1),min(ylim(2),CONFIG.ZOOM_TPOT)]);
set(gca,'XLimMode','auto')
hold on
grid on

xlabel(F.psal.long_name)
ylabel('Potential temperature')

% subplot(1,2,1)
% plot(F.temp.data(ideb:ifin,:)',F.pres.data(ideb:ifin,:)',Marker);
% hold on
% hl=plot(F.temp.data(ideb,:)',F.pres.data(ideb,:)',Marker);
% set(gca,'Ydir','reverse')
% hold on
% grid on


%keyboard
%nocycl=find(max(F.pres.data')==max(max(F.pres.data)))
nocycl=2
[F,Dim]=extract_profile_dim(F,Dim,'N_PROF',nocycl);

% supperpose la CTD de ref
if isempty(campaign_name)==0
    CAMPAIGN=load_configuration(['../VERIF_PROF1/config_campagne/' campaign_name ]);
end
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
if isempty(campaign_name)==0
    filename = fullfile( CAMPAIGN.DATA_DIRECTORY, strcat(CAMPAIGN.CAMPAGNE_MAT, CAMPAIGN.POSTFIX ));
    CTD_ALL = load( filename );
    
    %lat_campagne = CTD.LAT;
    %lon_campagne = CTD.LONG;
    dista1 = andoyer(CTD_ALL.LONG,CTD_ALL.LAT, F.longitude.data(:,:),F.latitude.data(:,:));
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
    CTD.latitude.data
    
    %---------------
    % Interpolation
    %---------------
    % default CONFIG
    PARAM.DIRECTION='D';
    PARAM.LIM_SURF_DEEP=500;
    PARAM.STEP_SURF=10;
    PARAM.STEP_DEEP=20;
    PARAM.DATATYPE='raw';
    PARAM.NB_PROF=1;
    PARAM.ZOOM=1000;
    PARAM.DATAREP='DIR_FTP';
    
    addpath('../VERIF_PROF1/util/')
    CTDi_on_pres = interpolation_m(CTD,F,'pres',PARAM); % interpolation of CTD variables on float presssure levels F.pres
    
    CTDi_on_ptmp0 = interpolation_m(CTD,F,'ptmp0',PARAM); % interpolation of CTD variables  on  float theta levels F.ptmp0
    
    hi=plot(CTDi_on_ptmp0.psal.data(:,:),CTDi_on_ptmp0.ptmp0.data(:,:),'.k');
end
%legend(hi,['CTD made at float launch ' F.platform_number.data(1,:)])
