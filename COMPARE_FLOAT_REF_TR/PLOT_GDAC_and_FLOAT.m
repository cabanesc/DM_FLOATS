% -========================================================
%   USAGE : PLOT_GDAC_and_FLOAT(floatname,dacname,varargin)
%   PURPOSE : find and plot all Argo profiles availble on GDAC FTP  and located close to a given Argo float (map on theta levels, time series on theta levels theta/S diagram)
% -----------------------------------
%   INPUT :
%    floatname  (char)  e.g. '690258'
%    dacname    (char) e.g.  'coriolis'
%
%   OPTIONNAL INPUT :
%    'DEPTH_MIN' (float)        1000 (default=1000)minimum depth allowed when selecting data on specified theta levels
%    'YEAR_MIN'  (float)        minimum year allowed for gdac data (default=1989)
%    'YEAR_MAX'  (float)        maximum year allowed for gdac data  (default=current year)
%    'ONTHETA'   (float)        =1 if plots are made on theta levels (TPOT_MIN-TPOT_MAX)
%                               =0 if plots are made on pressure levels (P_MIN P_MAX) (default =1)
%    'TPOT_MIN'  (float)        (default =3)
%    'TPOT_MAX'  (float)        (default =3.2)
%    'P_MIN'     (float)        (default =1800)
%    'P_MAX'     (float)        (default =2000)
%    'THE_BASES' (cell of chars) ex: {'CTD','ARGO'} (default) or {'CTD'} or {'ARGO'} :reference data taken into account
%    'CYCLE_MIN' (float)        minimum float's cycle to consider (default=1)
%    'CYCLE_MAX' (float)        maximum float's cycle to consider  (default=last cycle available)
%    'VEC_PSAL'  (1x2)          min and max PSAL values to scale colorbar for plots (e.g. [34.68 35.01])
%    'UPDATE'    (logical)      UPDATE=1 if the list of all Argo profiles availble
%                               on GDAC is updated each time this program is run, UPDATE =0 otherwise (default=1)
%
%    'DATATYPE'  (char)        'raw' if raw PARAM are plotted; 'adj' if PARAM_ADJ are plotted  (default is 'raw')
%    'PRINT'     (logical)      PRINT=1 is figure is saved  PRINT=0 otherwise
%                               (default=0)
% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  : created (2016) ccabanes
%
%   CALLED SUBROUTINES:
% ========================================================

function PLOT_GDAC_and_FLOAT(floatname,dacname,varargin)

close all
%init_path

CONFIG=load_configuration('config.txt');

CONFIG.pathwmobox=[CONFIG.DIR_OWC CONFIG.VERSION_OWC '/data/'];;

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

% To save list of nearby floats
CONFIG.file_nearby_floats=[CONFIG.plotpath '/list_' floatname '.mat'];

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
F=read_netcdf_allthefile([CONFIG.DIR_FTP dacname '/' floatname '/' floatname '_prof.nc'],NcVar);

% default CONFIG
CONFIG.DEPTH_MIN=1000;
CONFIG.YEAR_MIN=1989;
CONFIG.YEAR_MAX=str2num(datestr(date,'YYYY'));
CONFIG.ONTHETA=1;
CONFIG.TPOT_MIN=3.4;
CONFIG.TPOT_MAX=3.5;
CONFIG.P_MIN=1800;
CONFIG.P_MAX=2000;
CONFIG.CYCLE_MIN=1;
CONFIG.CYCLE_MAX=max(F.cycle_number.data);
CONFIG.THE_BASES={'CTD','ARGO'};
CONFIG.VEC_PSAL_AUTO='on';
CONFIG.UPDATE=1;
CONFIG.PRINT=0;
CONFIG.DATATYPE='raw';

% Input CONFIG
if isfield(s,'DEPTH_MIN')==1;CONFIG.DEPTH_MIN=s.DEPTH_MIN;end;
if isfield(s,'YEAR_MIN')==1;CONFIG.YEAR_MIN=s.YEAR_MIN;end;
if isfield(s,'YEAR_MAX')==1;CONFIG.YEAR_MAX=s.YEAR_MAX;end
if isfield(s,'ONETHETA')==1;CONFIG.ONTHETA=s.ONTHETA;end
if isfield(s,'TPOT_MIN')==1;CONFIG.TPOT_MIN=s.TPOT_MIN;end
if isfield(s,'TPOT_MAX')==1;CONFIG.TPOT_MAX=s.TPOT_MAX;end
if isfield(s,'P_MIN')==1;CONFIG.P_MIN=s.P_MIN;end
if isfield(s,'P_MAX')==1;CONFIG.P_MAX=s.P_MAX;end
if isfield(s,'THE_BASES')==1;CONFIG.THE_BASES=s.THE_BASES;end
if isfield(s,'CYCLE_MIN')==1;CONFIG.CYCLE_MIN=s.CYCLE_MIN;end
if isfield(s,'CYCLE_MAX')==1;CONFIG.CYCLE_MAX=s.CYCLE_MAX;end
if isfield(s,'UPDATE')==1;CONFIG.UPDATE=s.UPDATE;end
if isfield(s,'PRINT')==1;CONFIG.PRINT=s.PRINT;end
if isfield(s,'VEC_PSAL')==1;
    CONFIG.VEC_PSAL=s.VEC_PSAL;
    CONFIG.VEC_PSAL_AUTO='off';
end
if isfield(s,'DATATYPE')==1;CONFIG.DATATYPE=s.DATATYPE;end;


THE_BASES=CONFIG.THE_BASES;

% regional box

ideb=find(F.cycle_number.data==CONFIG.CYCLE_MIN);
ifin=find(F.cycle_number.data==CONFIG.CYCLE_MAX);

time_slot=[ideb(1):ifin(1)];


CONFIG.VEC_REG=[(floor(min(F.longitude.data(time_slot))*10)/10-10),ceil(max(F.longitude.data(time_slot))*10)/10,floor(min(F.latitude.data(time_slot))*10)/10,(ceil(max(F.latitude.data(time_slot))*10)/10+3)]

[topo,DimT]=read_netcdf_allthefile(CONFIG.FILE_TOPO);

if strcmp(CONFIG.DATATYPE,'adj') == 1
    thesrt='(adj)';
    isbest_value=2;
else
    isbest_value=0;
    thestr='(raw)';
end

% COMPARAISON AUX DONNEES ARGO TEMPS REEL.
%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,1,1)
title({[' REAL TIME and float data ' floatname ' ' thestr ' :'];[' salinity at theta levels ' num2str(CONFIG.TPOT_MIN) ' - ' num2str(CONFIG.TPOT_MAX)]});
contour_bathy(topo,CONFIG.VEC_REG)
hold on
box on
grid on

subplot(2,1,2)
hold on
box on
grid on
xlabel('year')
ylabel(['Salinity' ])

% trouve les flotteurs TR dans la zone
if ~exist(CONFIG.file_nearby_floats) || CONFIG.UPDATE==1
    index_file_name=[CONFIG.DIR_FTP_CORIOLIS '../ar_index_global_prof.txt'];
    index_data = read_index_prof(index_file_name);
    
    [liste_profile_files,liste_floats_files] = select_floats_from_index(index_data,'lonmin',CONFIG.VEC_REG(1), 'lonmax', CONFIG.VEC_REG(2), 'latmin',CONFIG.VEC_REG(3), 'latmax', CONFIG.VEC_REG(4) ,'ocean','A','exclude_type',{'852','851'},'datemin','20020101');
    save(CONFIG.file_nearby_floats,'liste_profile_files','liste_floats_files')
else
    load(CONFIG.file_nearby_floats);
    disp(['Number of  nearby floats selected : ' num2str(length(liste_floats_files))])
end

ctd_to_exclude=floatname;
[tabdac,r]=strtok(liste_floats_files,'/');
tabfloat=strtok(r,'/');
iio=ismember(tabfloat,ctd_to_exclude);
liste_floats_files(iio)=[];
CONFIG_GDAC=CONFIG;
CONFIG_GDAC.DIR_FTP=CONFIG_GDAC.DIR_FTP_CORIOLIS;
second_floatname='';
%for kfloat=1:length(liste_floats_files)
for kfloat=1:50
    disp(['nearby float:' num2str(kfloat) '/' num2str(length(liste_floats_files))]);
    subplot(2,1,1)
    plot_float_on_map(tabfloat{kfloat},tabdac{kfloat},topo,CONFIG_GDAC,'psal_mean','MarkerSize',20,'MarkerLine','off','Isbest',1);
    subplot(2,1,2)
    plot_float_on_time(tabfloat{kfloat},tabdac{kfloat},CONFIG_GDAC,'','MarkerFaceColor','b','MarkerSize',50,'MarkerLine','off','Isbest',1,'Number_floats',[kfloat,length(liste_floats_files)]);
    if strcmp(strtrim(tabfloat{kfloat}),strtrim(second_floatname))
        subplot(2,1,1)
        plot_float_on_map(tabfloat{kfloat},tabdac{kfloat},topo,CONFIG_GDAC,'psal_mean','MarkerSize',50,'MarkerLine','off','Isbest',1);
        subplot(2,1,2)
        plot_float_on_time(tabfloat{kfloat},tabdac{kfloat},CONFIG_GDAC,'','MarkerFaceColor','c','MarkerSize',50,'MarkerLine','off','Isbest',1,'Number_floats',[kfloat,length(liste_floats_files)]);
    end
end

thexlim=get(gca,'YLim');
subplot(2,1,2)
plot_float_on_time(floatname,dacname,CONFIG,'','MarkerFaceColor','b','MarkerSize',50,'Isbest',isbest_value)
plot_float_on_time(floatname,dacname,CONFIG,'','MarkerSize',50,'Isbest',isbest_value)


subplot(2,1,1)
plot_float_on_map(floatname,dacname,topo,CONFIG,'psal_mean','MarkerSize',50,'Isbest',isbest_value)
if strcmp(CONFIG.VEC_PSAL_AUTO,'on')==0
    caxis(CONFIG.VEC_PSAL)
else
    caxis(thexlim)
end
colorbar
titplt=[floatname '_and_CURR_at_theta.png'];
if CONFIG.PRINT==1
    eval(['print -dpng ' CONFIG.plotpath  titplt]);
end

