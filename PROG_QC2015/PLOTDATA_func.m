% -========================================================
%   USAGE :   PLOTDATA_FUNC(floatname,dacname,varargin)
%   PURPOSE : plot data from a given float
% -----------------------------------
%   INPUT :
%    floatname  (char)  e.g. '690258'
%    dacname    (char) e.g.  'coriolis'
%
%   OPTIONNAL INPUT :
%    'DATATYPE' (char)  'raw' if raw PARAM are plotted 'adj' if PARAM_ADJ are plotted  (default is 'raw')
%    'DATAREP'   (char)   default 'DIR_FTP'  directory where nc files are looked for  or 'DIR_DM_FILES'
%    'USEFLAG'  (logical) ==1 if flag are taken into account for plot ==0
%    if not (default 0)
%    'PRINT'     (logical)      PRINT=1 is figure is saved  PRINT=0 otherwise
%                               (default=1)
%    'MAKEPAUSE'  (logical)     MAKEPAUSE=1 : PAUSE at each profile when
%    plotting the theta/s diagram (default=0)
%    'PLOT_INVDENS' (logical)   PLOT_INVDENS=1, create and display all the figures of density inversion (default)
%                                           =0, do not create  the figures of density inversion
% -----------------------------------
%   OUTPUT :
% -----------------------------------
%   HISTORY  :  created: V. Thierry - N. David - Juillet 2007
%               revised: C. Lagadec
%               revised: ccabanes - 2016
%   CALLED SUBROUTINES:
% ========================================================

%function PLOTDATA_FUNC(float,dacname,datatype,inuseflag,park_press,prof_press)
function PLOTDATA_FUNC(floatname,dacname,rdir,varargin)
close all

%init_path % ==> Commented by T. Reynaud 07/09/2020

CONFIG=load_configuration([rdir,'config.txt']);


n=length(varargin);

if n/2~=floor(n/2)
    error('check the imput arguments')
end

f=varargin(1:2:end);
c=varargin(2:2:end);
s = cell2struct(c,f,2);

% default CONFIG
PARAM.DATATYPE='raw';
PARAM.USEFLAG='n';
PARAM.PRINT=1;
PARAM.DATAREP='DIR_FTP';
PARAM.MAKEPAUSE=0;
PARAM.PLOT_INVDENS=1;

% Input CONFIG
if isfield(s,'DATATYPE')==1;PARAM.DATATYPE=s.DATATYPE;end;
if isfield(s,'USEFLAG')==1;PARAM.USEFLAG=s.USEFLAG;end;
if isfield(s,'PRINT')==1;PARAM.PRINT=s.PRINT;end;
if isfield(s,'DATAREP')==1;PARAM.DATAREP=s.DATAREP;end;
if isfield(s,'MAKEPAUSE')==1;PARAM.MAKEPAUSE=s.MAKEPAUSE;end;
if isfield(s,'PLOT_INVDENS')==1;PARAM.PLOT_INVDENS=s.PLOT_INVDENS;end;
float=str2num(floatname);

%------------------------------------------------------
% lecture des meta pour recuperer park_press et prof_press
%------------------------------------------- -----------

M = read_netcdf_allthefile([CONFIG.DIR_FTP dacname '/' floatname '/' floatname '_meta.nc']);
[isfound]=findstr_tab(M.config_parameter_name.data,'ParkPressure_dbar');
if sum(isfound)==1
    park_press=median(M.config_parameter_value.data(:,isfound));
    
else
    park_press=1000;
end

disp(['Parking presssure: ' num2str(park_press)])
[isfound]=findstr_tab(M.config_parameter_name.data,'ProfilePressure_dbar');
if sum(isfound)==1
    prof_press=median(M.config_parameter_value.data(:,isfound));
else
    prof_press=2000;
end
disp(['Profile presssure: ' num2str(prof_press)])

if strcmp(PARAM.USEFLAG,'n') == 1
    useflag=0;
    titflag='flagnotused';
else
    useflag=1;
    titflag='flagused';
end

if strcmp(PARAM.DATATYPE,'adj') ==1 & strcmp(PARAM.DATAREP,'DIR_DM_FILES')==1
    titdata='DM';
    titsavedata='';
    FLOAT_SOURCE_NETCDF = CONFIG.(PARAM.DATAREP);
    
    %FLOAT_SOURCE_NETCDF = CONFIG.DIR_OUT;
    %syst_config = load_configuration( 'config_all_dm.txt' );
    %dacname='';
    titflag=[titflag '_dm'];
    datatype=2;
elseif strcmp(PARAM.DATATYPE,'adj') == 1 & strcmp(PARAM.DATAREP,'DIR_FTP')==1
    datatype=2;
    titdata='RT';
    titsavedata='';
    %FLOAT_SOURCE_NETCDF = CONFIG.DIR_FTP;
    FLOAT_SOURCE_NETCDF = CONFIG.(PARAM.DATAREP);
    %syst_config = load_configuration( 'config_deep.txt' );
elseif strcmp(PARAM.DATATYPE,'raw')==1
    datatype=1;
    titdata='RT';
    titsavedata='';
    %FLOAT_SOURCE_NETCDF = CONFIG.DIR_FTP;
    FLOAT_SOURCE_NETCDF = CONFIG.(PARAM.DATAREP);
    %syst_config = load_configuration( 'config_deep.txt' );
end



floatnum=floatname;

display(['  Float : ' floatnum]);


repnc=[FLOAT_SOURCE_NETCDF '/' dacname '/' floatnum '/profiles']

plotpath = [CONFIG.DIR_PLOT '/preliminaire/' floatnum '/'];
if ~exist(plotpath,'dir')
    mkdir(plotpath);
end

% lecture des fichiers TECH pour recuperer certaines variables
FILENAME_TECH = [CONFIG.DIR_FTP dacname '/' floatname '/' floatname '_tech.nc'];
if exist(FILENAME_TECH)
    T = read_netcdf_allthefile(FILENAME_TECH);
    % trouve les pressions de surface
    isurf=find(findstr_tab(cellstr(T.technical_parameter_name.data),'PRES_SurfaceOffset')==1);
    pres_name = unique(cellstr(T.technical_parameter_name.data(isurf,:)));
    ivolt=find(findstr_tab(cellstr(T.technical_parameter_name.data),'VOLTAGE_BatteryPumpStartProfile_volts')==1);
    batt_name = unique(cellstr(T.technical_parameter_name.data(ivolt,:)));
    ipump=find(findstr_tab(cellstr(T.technical_parameter_name.data),'NUMBER_PumpActionsDuringAscentToSurface_COUNT')==1);
    pump_name = unique(cellstr(T.technical_parameter_name.data(ipump,:)));
    ivalve=find(findstr_tab(cellstr(T.technical_parameter_name.data),'NUMBER_ValveActionsDuringDescentToPark_COUNT')==1);
    valve_name = unique(cellstr(T.technical_parameter_name.data(ivalve,:)));
    ivalve=find(findstr_tab(cellstr(T.technical_parameter_name.data),'NUMBER_ValveActionsDuringDescentToPark_COUNT')==1);
    valve_name = unique(cellstr(T.technical_parameter_name.data(ivalve,:)));
    
    cycl=T.cycle_number.data(isurf);
    cyc1=(cycl==1);
    
    figure
    
    isplot=0;
    if isempty(pres_name)==0&length(pres_name)==1
        subplot(3,1,1)
        
        plot(T.cycle_number.data(isurf),str2num(T.technical_parameter_value.data(isurf,:)),'+-')
        hold on
        box on
        grid on
        xlabel('cycle_number','interpreter','none')
        ylabel('dbar')
        title(pres_name{1},'interpreter','none')
        isplot=1;
    end
    if isempty(batt_name)==0&length(batt_name)==1
        subplot(3,1,2)
        
        plot(T.cycle_number.data(ivolt),str2num(T.technical_parameter_value.data(ivolt,:)),'+-')
        hold on
        box on
        grid on
        xlabel('cycle_number','interpreter','none')
        ylabel('volt')
        title(batt_name{1},'interpreter','none')
        isplot=1;
        
    end
    if isempty(pump_name)==0&isempty(valve_name)==0
        subplot(3,1,3)
        
        plot(T.cycle_number.data(ipump),str2num(T.technical_parameter_value.data(ipump,:)),'+-')
        hold on
        plot(T.cycle_number.data(ivalve),str2num(T.technical_parameter_value.data(ivalve,:)),'o-r')
        hold on
        box on
        grid on
        xlabel('cycle_number','interpreter','none')
        ylabel('COUNT')
        title({'NUMBER_PumpActionsDuringAscentToSurface_COUNT (+)'; 'NUMBER_ValveActionsDuringDescentToPark_COUNT (o)'},'interpreter','none')
        isplot=1;
    end
    if isplot==1
        fileout=[CONFIG.DIR_PLOT 'preliminaire/' floatname '/' floatname '_surface_pres.pdf'];
        if exist([CONFIG.DIR_PLOT 'preliminaire/' floatname '/'])==0
            mkdir([CONFIG.DIR_PLOT 'preliminaire/' floatname '/'])
        end
        %keyboard
        set(gcf,'Position',[1 46 660 800],'paperPositionMode','auto');
        print(fileout,'-dpdf')
    end
end
%------------------------------------------------------
% lecture de tous les monoprofifs du flotteur
% restitution dans structure FLm 
%------------------------------------------- -----------
[file_list]=select_float_files_on_ftp(floatname,dacname,FLOAT_SOURCE_NETCDF,'C',0);
vertical_sampling_scheme='Primary sampling';
Param='';
[FLm,DimL,file_list]=create_multi_from_filelist(floatname,dacname,FLOAT_SOURCE_NETCDF,file_list,vertical_sampling_scheme,Param);
FLm=replace_fill_bynan(FLm); % add 23/01/2024


% lecture des fichiers TECH aux pour recuperer certaines variables
%------------------------------------------------------
isplot_ice=0;
FILENAME_TECH_AUX = [CONFIG.DIR_FTP dacname '/' floatname '/auxiliary/' floatname '_tech_aux.nc'];
if exist(FILENAME_TECH_AUX)
    T = read_netcdf_allthefile(FILENAME_TECH_AUX);
    iice=find(findstr_tab(cellstr(T.technical_parameter_name.data),'TECH_FLAG_IceDetection_NUMBER')==1);
    icedetect_name = unique(cellstr(T.technical_parameter_name.data(iice,:)));
    iiceact=find(findstr_tab(cellstr(T.technical_parameter_name.data),'TECH_FLAG_IceAlgorithmActivated_LOGICAL')==1);
    iceactive_name = unique(cellstr(T.technical_parameter_name.data(iiceact,:)));
    itrans=find(findstr_tab(cellstr(T.technical_parameter_name.data),'TECH_FLAG_TransmissionDelayed_NUMBER')==1);
    trans_name=unique(cellstr(T.technical_parameter_name.data(itrans,:)));
    
    full_cycle=[1:max(T.cycle_number.data(iice))];
    empty_cycle=~ismember(full_cycle,T.cycle_number.data(iice));
    empty_cycle_num=find(~ismember(full_cycle,T.cycle_number.data(iice)));
   
    m=1;
    if isempty(iceactive_name)==0
        if sum(str2num(T.technical_parameter_value.data(iiceact,:)))>0 % si au moins un cycle a l'algo ice active
            figure
            t1 = tiledlayout(3,1,'TileSpacing','compact');
            set(gcf,'Position',[1 46 660 800],'paperPositionMode','auto');
            nexttile
            %subplot(3,1,1)
            hold on
            b=bar(full_cycle,empty_cycle*4);
            b.FaceColor=[0.9 0.9 0.9];
            b.LineStyle='none';
            %b.EdgeColor=[1 1 1];
            box on
            grid on
            xlabel('cycle_number','interpreter','none')
            ylabel('Number (logical)')
            xl=get(gca,'XLim');
            set(gca,'XLim',[xl(1),max(FLm.cycle_number.data)]);
            if isempty(trans_name)==0&length(trans_name)==1
                %subplot(3,1,1)
                %nexttile
                trans_num=str2num(T.technical_parameter_value.data(itrans,:));
                trans_num(ismember(T.cycle_number.data(itrans),empty_cycle_num))=NaN;
                if isempty(trans_num)==0
                    p=plot(T.cycle_number.data(itrans),trans_num,'.-');
                    p.MarkerSize=8;
                    m=max(m,max(trans_num));
                    isplot_ice=isplot_ice+1;
                end
            end
            if isempty(icedetect_name)==0&length(icedetect_name)==1
                ice_num=str2num(T.technical_parameter_value.data(iice,:));
                if isempty(ice_num)==0
                    p1=plot(T.cycle_number.data(iice),str2num(T.technical_parameter_value.data(iice,:)),'o-');
                    %p1.Color=[0.301 0.745 0.933];
                    p1.Color=[0 0.4470 0.7410];
                    m=max(m,max(ice_num));
                    %str=(icedetect_name{1},'interpreter','none')
                    isplot_ice=isplot_ice+2;
                end
            end
            set(gca,'YLim',[0 m])
            if isplot_ice>0
                if isplot_ice==2
                    hl=legend([p1],'Ice detection number : (0:no detection, 1:ISA , 2:sat mask, 3:ISA+sat mask, 4:ascent hanging )','location','southoutside');

                    title([floatname ': Ice detection'],'interpreter','none');
                elseif isplot_ice==1
                    hl=legend([p],'Transmission delayed number : (0:no delay, 1:all data delayed, 2:some data delayed)','location','southoutside');

                    title([floatname ': Transmission delayed number'],'interpreter','none');
                elseif isplot_ice==3
                    hl=legend([p,p1],{'Transmission delayed number : (0:no delay, 1:all data delayed, 2:some data delayed)';'Ice detection number : (0:no detection, 1:ISA , 2:sat mask, 3:ISA+sat mask, 4:ascent hanging '},'location','southoutside');
                    title([floatname ': Ice detection and related parameters'],'interpreter','none');
                end
                
            end
        end
    end
end
%------------------------------------------------------
% Profils du flotteur: mise dans variables
%------------------------------------------- -----------


read_ARGO_tous_mono3_1

[nlev,ncyc]=size(temp);

if strcmp(pltoxy ,'n')
    doxy_qc_num(1:nlev,1:ncyc) = ' ';
    doxy(1:nlev,1:ncyc)    = NaN;
    psat(1:nlev,1:ncyc)    = NaN;
end

if  useflag == 1
    inan=find(psal_qc_num==3 | psal_qc_num==4 | pres_qc_num == 3 | pres_qc_num == 4);
    psal(inan)= NaN;
    sig0(inan)= NaN;
    
    inan=find(temp_qc_num==3 | temp_qc_num == 4 | pres_qc_num == 3 | pres_qc_num == 4);
    temp(inan) = NaN;
    sig0(inan) = NaN;
    tpot(inan) = NaN;
    
    inan=find(pres_qc_num==3 | pres_qc_num==4);
    pres(inan)=NaN;
    % si P=2047 P est flaguee a 4, pas T et S
    temp(inan) = NaN;
    psal(inan) = NaN;
    tpot(inan) = NaN;
    sig0(inan) = NaN;
    
    inan=find(doxy_qc_num==3 | doxy_qc_num==4 | pres_qc_num == 3 | pres_qc_num == 4);
    doxy(inan)=NaN;
    doxy(doxy>370)=NaN;
end

pres = pres';
temp = temp';
psal = psal';
doxy = doxy';
sig0 = sig0';
tpot = tpot';

map=jet(ncyc);

% calcul TPOT(ref a 1000 dbar)
% tpot et sig0 sont calcules dans read_ARGO_mono3.1
% calcul SIG1
% calcul Saturation O2

tpot1 = sw_ptmp(psal,temp,pres,1000);
sig1  = sw_pden(psal,temp,pres,1000)-1000;

if  strcmp(pltoxy,'o')
    sato2 = sw_satO2(psal,temp);
    [sato2b]=convert_oxygen(sato2,'mL/L','mumol/kg',sig0);
    psat=100*doxy./sato2b;
end

% lecture et trace des parametres techniques du fichier .tech
% (batterie et pression de surface)
% -----------------------------------------------------------

%    trait_param_tech(dacname,floatnum,plotpath,cycnum)

% trace des donnees brutes en fonction de la pression
% TEMP,PSAL,SIG0,TPOT
% PSAT et DOXY si oxygene
% ---------------------------------------------------
if strcmp(PARAM.DATATYPE,'adj')
    plot_dbrut(plotpath,floatnum,psal,pres,'PSAL_ADJ','PRES','PSAL (PSU)','PRES (dbar)',cycnum,lon,titflag)
    plot_dbrut(plotpath,floatnum,temp,pres,'TEMP_ADJ','PRES','TEMP (deg.celsius)','PRES (dbar)',cycnum,lon,titflag)
    plot_dbrut(plotpath,floatnum,tpot,pres,'TPOT_ADJ','PRES','TPOT (deg.celsius)','PRES (dbar)',cycnum,lon,titflag)
    plot_dbrut(plotpath,floatnum,sig0,pres,'SIG0_ADJ','PRES','SIG0 ','PRES (dbar)',cycnum,lon,titflag)
    if  strcmp(pltoxy,'o')
        plot_dbrut(plotpath,floatnum,doxy,pres,'DOXY_ADJ','PRES','DOXY (mumol/kg)','PRES (dbar)',cycnum,lon,titflag)
        plot_dbrut(plotpath,floatnum,psat,pres,'O2sat_ADJ','PRES','O2 % saturation','PRES (dbar)',cycnum,lon,titflag)
    end
else
    plot_dbrut(plotpath,floatnum,psal,pres,'PSAL','PRES','PSAL (PSU)','PRES (dbar)',cycnum,lon,titflag)
    plot_dbrut(plotpath,floatnum,temp,pres,'TEMP','PRES','TEMP (deg.celsius)','PRES (dbar)',cycnum,lon,titflag)
    plot_dbrut(plotpath,floatnum,tpot,pres,'TPOT','PRES','TPOT (deg.celsius)','PRES (dbar)',cycnum,lon,titflag)
    plot_dbrut(plotpath,floatnum,sig0,pres,'SIG0','PRES','SIG0 ','PRES (dbar)',cycnum,lon,titflag)
    if  strcmp(pltoxy,'o')
        plot_dbrut(plotpath,floatnum,doxy,pres,'DOXY','PRES','DOXY (mumol/kg)','PRES (dbar)',cycnum,lon,titflag)
        plot_dbrut(plotpath,floatnum,psat,pres,'O2sat','PRES','O2 % saturation','PRES (dbar)',cycnum,lon,titflag)
    end
end

% trace des niveaux verticaux en fonction de la pression
% (mise sur des niveaux identiques)
% -----------------------------------------------------------------------
%plante  si lance en VPN
% pres_qc_num = pres_qc_num';
% if length(cycnum)>2
% plot_levelvert(plotpath,floatnum,cycnum,pres,pres_qc_num,titflag,titsavedata)
% end
disp('ok')

% trace des inversions de densite (test 14 de la doc ARGO)
% --------------------------------------------------------

if PARAM.PLOT_INVDENS==1
    plot_inverdens(CONFIG,floatnum,cycnum,pres,temp,psal,titflag)
end

% Diagramme T/S, theta/S, theta/doxy
%-----------------------------------

h2=figure;
orient landscape
set(gca,'Fontsize',16);
hold on
for iprof=1:ncyc
    plot(psal(iprof,:),temp(iprof,:),'color',map(iprof,:),'marker','.','linewidth',2);
end
vxmax=max(ceil(psal(:)*100))/100;
vxmin=min(floor(psal(:)*100))/100;
vymax=max(ceil(temp(:)*100))/100;
vymin=min(floor(temp(:)*100))/100;
if ~isnan(vxmax*vxmin*vymax*vymin)
    set(gca,'xlim',[vxmin vxmax],'ylim',[vymin vymax]);
end
if strcmp(PARAM.DATATYPE,'adj')
    ylabel('Temperature Adj.')
    xlabel('Salinity Adj.')
else
    ylabel('Temperature')
    xlabel('Salinity')
end
grid on
box on
title([floatnum ', T/S diagram']);
if PARAM.PRINT==1
    set(gcf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]);
    %eval(['print -depsc2 ' plotpath floatnum '_TS_raw_' titflag titsavedata '.eps']);
    eval(['print -f' num2str(h2.Number) ' -dpng ' plotpath floatnum '_TS_raw_' titflag titsavedata '.png']);
end

%------
h3=figure;
orient landscape
set(gca,'Fontsize',16);
hold on
for iprof=1:ncyc
    plot(psal(iprof,:),tpot(iprof,:),'color',map(iprof,:),'marker','.','linewidth',2);
end
if ~isnan(min(tpot(:)).*max(tpot(:)))
    ctpot=[floor(min(tpot(:)))-1:0.1:ceil(max(tpot(:)))+1]';
    cpsal=[floor(10*min(psal(:)))/10-0.1:0.1:ceil(max(psal(:))*10)/10+0.1]';
    [tabct,tabcp]=meshgrid(ctpot,cpsal);
    [~,csig]=swstat90(tabcp,tabct,0);
    if size(csig,1)~=size(tabcp,1)
        csig=csig';
    end
    [c,h]=contour(tabcp,tabct,csig,[20:0.5:35],'k');
    %keyboard
    clabel(c,h,'LabelSpacing',144*2,'FontSize',10,'FontWeight','normal')
    if strcmp(PARAM.DATATYPE,'adj')
        ylabel('Potential Temp. Adj. (ref. to 0db)')
        xlabel('Salinity Adj.')
    else
        ylabel('Potential Temp. (ref. to 0db)')
        xlabel('Salinity')
    end
    vxmax=max(ceil(psal(:)*100))/100;
    vxmin=min(floor(psal(:)*100))/100;
    vymax=max(ceil(tpot(:)*100))/100;
    vymin=min(floor(tpot(:)*100))/100;
    if ~isnan(vxmax*vxmin*vymax*vymin)
        set(gca,'xlim',[vxmin vxmax],'ylim',[vymin vymax]);
    end
end
grid on
box on
title([floatnum ', theta/S diagram']);

if PARAM.PRINT==1
    set(gcf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]);
    %eval(['print -depsc2 ' plotpath floatnum '_thetaS_' titflag titsavedata '.eps']);
    eval(['print -f' num2str(h3.Number) ' -dpng ' plotpath floatnum '_thetaS_' titflag titsavedata '.png']);
end

%------
h4=figure;
%keyboard
orient landscape
set(gca,'Fontsize',16);
hold on
pres_tab=pres;
psal_tab=psal;
tpot_tab=tpot;
psal_tab(pres_tab<0)=NaN;
tpot_tab(pres_tab<0)=NaN;
for iprof=1:ncyc
    plot(psal_tab(iprof,:),tpot_tab(iprof,:),'color',map(iprof,:),'marker','.','linewidth',2);
    if(PARAM.MAKEPAUSE==1)
        disp(num2str(cycnum(iprof)));
        % pause
    end
end
psal_tab(pres_tab<1500)=NaN;
tpot_tab(pres_tab<1500)=NaN;
ctpot=[floor(min(tpot_tab(:)))-1:0.05:ceil(max(tpot_tab(:)))+1]';
cpsal=[floor(10*min(psal_tab(:)))/10-0.1:0.05:ceil(max(psal_tab(:))*10)/10+0.1]';
[tabct,tabcp]=meshgrid(ctpot,cpsal);
[~,csig]=swstat90(tabcp,tabct,0);
if size(csig,1)~=size(tabcp,1)
    csig=csig';
end
%[c,h]=contour(tabcp,tabct,csig,[25:0.05:35],'k');
%keyboard
%clabel(c,h,'LabelSpacing',144,'FontSize',10,'FontWeight','demi')
%[c,h]=contour(tabcp,tabct,csig,[25:0.01:35],'k:');
if strcmp(PARAM.DATATYPE,'adj')
    ylabel('Potential Temp. Adj. (ref. to 0db)')
    xlabel('Salinity Adj.')
else
    ylabel('Potential Temp. (ref. to 0db)')
    xlabel('Salinity')
end
%keyboard
vxmax=max(ceil(psal_tab(:)*100))/100;
vxmin=min(floor(psal_tab(:)*100))/100;
vymax=max(ceil(tpot_tab(:)*100))/100;
vymin=min(floor(tpot_tab(:)*100))/100;
if ~isnan(vxmax*vxmin*vymax*vymin)
    set(gca,'xlim',[vxmin vxmax],'ylim',[vymin vymax]);
    
    %set(gca,'xlim',[cpsal(1) cpsal(end)],'ylim',[ctpot(1) ctpot(end)]);
    if max(max(csig)-min(min(csig)))<4
        [c,h]=contour(tabcp,tabct,csig,[0:0.05:50],'k');
        clabel(c,h,'LabelSpacing',2*144,'FontSize',10,'FontWeight','normal')
    else
        [c,h]=contour(tabcp,tabct,csig,[0:0.5:50],'k');
        clabel(c,h,'LabelSpacing',144*2,'FontSize',10,'FontWeight','normal')
    end
    if max(max(csig)-min(min(csig)))<2
        [c,h]=contour(tabcp,tabct,csig,[0:0.01:50],'k:');
    end
end
grid on
box on
title([ floatnum ', theta/S diagram (P>1500db)']);
if PARAM.PRINT==1
    set(gcf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]);
    %eval(['print -depsc2 ' plotpath floatnum '_thetaS_zoom_' titflag titsavedata '.eps']);
    eval(['print -f' num2str(h4.Number) ' -dpng ' plotpath floatnum '_thetaS_zoom_' titflag titsavedata '.png']);
end
%keyboard

if strcmp(pltoxy,'o')
    h5=figure;
    orient landscape
    set(gca,'Fontsize',18);
    hold on
    
    for iprof=1:ncyc
        plot(doxy(iprof,:),tpot(iprof,:),'color',map(iprof,:),'marker','.','linewidth',2);
    end
    ctpot=[floor(min(tpot(:)))-1:0.1:ceil(max(tpot(:)))+1]';
    cdoxy=[min(doxy(:))-10:10:max(doxy(:))+10]';
    [tabct,tabcp]=meshgrid(ctpot,cdoxy);
    if strcmp(PARAM.DATATYPE,'adj')
        ylabel('Potential Temp. Adj.(ref. to 0db)')
        xlabel('Oxygen Adj.')
    else
        ylabel('Potential Temp. (ref. to 0db)')
        xlabel('Oxygen')
    end
    vxmax=max(ceil(doxy(:)*100))/100;
    vxmin=min(floor(doxy(:)*100))/100;
    vymax=max(ceil(tpot1(:)*100))/100;
    vymin=min(floor(tpot1(:)*100))/100;
    if ~isnan(vxmax*vxmin*vymax*vymin)
        set(gca,'xlim',[vxmin vxmax],'ylim',[vymin vymax]);
    end
    %set(gca,'xlim',[cdoxy(1) cdoxy(end)],'ylim',[ctpot(1) ctpot(end)]);
    grid
    title(['Float WMO ' floatnum ', theta/O2 diagram']);
    if PARAM.PRINT==1
        set(gcf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]);
        eval(['print -depsc2 ' plotpath floatnum '_thetaoxy_' titflag titsavedata '.eps']);
        eval(['print -dpng ' plotpath floatnum '_thetaoxy_' titflag titsavedata '.png']);
    end
end

if strcmp(pltoxy,'o')
    tabparam = ['PSAL';'SIG0';'TPOT';'DOXY';'PSAT'];
    tpas      = [0.05,0.1,0.5,10,10];
    tcoef     = [10,10,1,0.1,0.1];
else
    tabparam = ['PSAL';'SIG0';'TPOT';'PSAT'];
    tpas      = [0.05,0.1,0.5,10];
    tcoef     = [10,10,1,0.1];
end

if strcmp(pltoxy,'o')
    ncas=5;
else
    ncas = 3;
end

for icas=1:ncas
    
    pas = tpas (icas);
    coef = tcoef(icas);
    if icas == 1 %PSAL
        tabval=psal;
        nomval=tabparam(icas,:);
        valini=psal;
        nomvalini='PSAL';
    elseif icas == 2 %SIG0
        tabval=sig0;
        nomval=tabparam(icas,:);
    elseif icas == 3 %TPOT
        tabval=tpot;
        nomval=tabparam(icas,:);
        valini=temp;
        nomvalini='TEMP';
    elseif icas == 4 %DOXY
        tabval=doxy;
        nomval=tabparam(icas,:);
        valini=doxy;
        nomvalini='DOXY';
    elseif icas == 5 %PSAT
        tabval=psat;
        nomval=tabparam(icas,:);
    end
    
    vmin=min(floor(tabval(:)*coef)/coef);
    vmax=max(ceil(tabval(:)*coef)/coef);
    tabcont=(vmin:pas:vmax);
    
    % Profondeur couche de melange
    clear tabmld
    critere='dsig&dt';critval=[0.03 0.2];
    
    for iprf=1:length(cycnum)
        [mlde,mlpt,mlps,mlpd]=calmld_val(sig0(iprf,:),pres(iprf,:),tpot(iprf,:),psal(iprf,:),critere,critval);
        tabmld(iprf)=mlde;
    end
    
    %keyboard
    tabtime=(juld-juld(1))*ones(1,size(pres,2));
    if length(cycnum)>2
        if (sum(sum(isnan(sig0)))<numel(sig0))&(sum(sum(isnan(tabval)))<numel(tabval))
            [hf] = pcolor_argodata_sig(cycnum,tabtime(:,1),pres',sig0',tabval',nomval,'interp');
        else
            hf=figure;
            grid on
            box on
        end
        %    plot(tabtime,tabmld','w-','linewidth',2);
        %    plot(cycnum,tabmld','w-','linewidth',2);
        
        orient landscape
        set(gca,'fontsize',18)
        %set(gca,'ydir','reverse');
        %colorbar
        %caxis([vmin vmax]);
        vtick=get(gca,'xtick');
        [ucy,ia]=unique(FLm.cycle_number.data);
        ujuld=juld(ia);
        datetic=interp1(ucy,ujuld,vtick);
        vticklabel=strcat(datestr(datenum(1950,1,1)+double(datetic),12));
        
        %vticklabel=datestr(datenum(1950,1,1)+vtick+juld(1),12);
        set(gca,'xtick',vtick,'xticklabel',vticklabel);
        
        set(gca,'xtick',vtick,'xticklabel',vticklabel);
        if strcmp(PARAM.DATATYPE,'adj')
            title([ floatnum ' - ' nomval ' ADJUSTED']);
        else
            title([ floatnum ' - ' nomval]);
        end
        ylabel('Pressure (db)')
        xlabel('Date')
        if PARAM.PRINT==1
            set(gcf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]);
            %['print -f' num2str(hf.Number) ' -dpng ' plotpath floatnum '_' nomval '_interp_' titflag titsavedata '.png']
            eval(['print -f' num2str(hf.Number) ' -dpng ' plotpath floatnum '_' nomval '_interp_' titflag titsavedata '.png'])
            %eval(['print -depsc2 ' plotpath floatnum '_' nomval '_interp_' titflag titsavedata '.eps'])
        end
        if icas ~= 2 && icas ~= 5
            
            if (sum(sum(isnan(valini)))<numel(valini))
                [hf] = pcolor_argodata(cycnum,pres',valini',nomvalini,'flat');
            else
                hf=figure;
                grid on
                box on
            end
            ylabel('Pressure (db)')
            xlabel('Cycle number')
            if strcmp(PARAM.DATATYPE,'adj')
                title([ floatnum ' - ' nomvalini ' ADJUSTED']);
                
                
            else
                title([ floatnum ' - ' nomvalini]);
                
            end
            %eval(['print -depsc2 ' plotpath floatnum '_' nomvalini '_raw_' titflag titsavedata '.eps']);
            %eval(['print -dpng ' plotpath floatnum '_' nomvalini '_raw_' titflag titsavedata '.png']);
        end
    end
end

%end
if isplot_ice>1 & length(cycnum)>2
    figure(2)
    %subplot(3,1,2)
    nexttile
    pres_surf=pres;
    temp_surf=temp;
    temp_surf(pres_surf>200)=NaN;
    pcolor_argodata_ice(cycnum,pres_surf',temp_surf','TEMP','flat');
    set(gca,'YLim',[0 50])
    title('Temperature measured in the upper 50-db layer')
    ylabel('Pressure(db)')
end
% carte 1 (positions)
% -------------------
%keyboard
zone_visu=[floor(min(lat))-0.2 ceil(max(lat))+0.2  floor(min(lon))-0.2 ceil(max(lon))+0.2];
reso='LR';proj='miller';
[hf,ha]=fct_pltmap_traj(zone_visu,reso,proj,park_press,prof_press);

%keyboard
for ii=2:length(lon)-1
    m_plot(lon(ii),lat(ii),'color',map(ii,:),'marker','o','markerfacecolor',map(ii,:),'markersize',7)
end
lonArrow = lon;
latArrow= lat;
tailleFleche = 0.761 * 0.005 * (zone_visu(4) -zone_visu(3))/16;
[xArrow, yArrow] = m_ll2xy(lonArrow, latArrow);
FLmc=format_flags_char2num(FLm);
for ii=1:length(lon)-1
    if FLmc.position_qc.data(ii)==8
        plot_arrow2(xArrow(ii), yArrow(ii), xArrow(ii+1), yArrow(ii+1), ...
            tailleFleche, tailleFleche, 'linestyle', '--','Color',[0.4 0.4 0.4],'FaceColor',[1 1 1],'EdgeColor',[0.4 0.4 0.4]);           
    elseif FLmc.position_qc.data(ii)==4|FLmc.position_qc.data(ii)==3
         plot_arrow2(xArrow(ii), yArrow(ii), xArrow(ii+1), yArrow(ii+1), ...
            tailleFleche, tailleFleche, 'linestyle', '-','Color',[0.4 0.4 0.4],'FaceColor',[1 0 0],'EdgeColor',[0.4 0.4 0.4]);
    elseif FLmc.position_qc.data(ii)==1|FLmc.position_qc.data(ii)==2
         plot_arrow2(xArrow(ii), yArrow(ii), xArrow(ii+1), yArrow(ii+1), ...
            tailleFleche, tailleFleche, 'linestyle', '-','Color',[0.4 0.4 0.4],'FaceColor',[0.4 0.4 0.4],'EdgeColor',[0.4 0.4 0.4]);
    else
        arrowHdl4 = plot_arrow2(xArrow(ii), yArrow(ii), xArrow(ii+1), yArrow(ii+1), ...
        0, 0, 'linestyle', '--','Color',[0.4 0.4 0.4],'FaceColor',[0.4 0.4 0.4],'EdgeColor',[0.4 0.4 0.4]);
    end
end


for ii=[5:5:length(lon)]
    j=m_text(lon(ii),lat(ii),int2str(cycnum(ii)));
    set(j,'color',map(ii,:),'fontsize',10,'hor','cen','VerticalAlignment','cap');
end
%m_plot(lon,lat,'k-')
ii=1;
hl1=m_plot(lon(ii),lat(ii),'color',map(ii,:),'marker','d','markerfacecolor',map(ii,:),'markersize',8);

ii=length(lon);
hl2=m_plot(lon(ii),lat(ii),'color',map(ii,:),'marker','s','markerfacecolor',map(ii,:),'markersize',8);

%htl=legend(hl,'First prof.','Last prof.','bestinside');
% pour legende uniquement
arrowHdl1= m_plot(-600,-600,'linestyle', '--','color',[0.4 0.4 0.4],'marker','>','markerfacecolor',[1 1 1],'markersize',5);
arrowHdl2= m_plot(-600,-600,'linestyle', '-','color',[0.4 0.4 0.4],'marker','>','markerfacecolor',[1 0 1],'markersize',5);
arrowHdl3= m_plot(-600,-600,'linestyle', '-','color',[0.4 0.4 0.4],'marker','>','markerfacecolor',[0.4 0.4 0.4],'markersize',5);
htl=legend([hl1 hl2 arrowHdl3 arrowHdl2 arrowHdl1],{'First prof.','Last prof.','QC=1,2','QC=3,4','QC=8'},'Location','northeast','Orientation','vertical');
%htl=legend([hl1 hl2 ],{'First prof.','Last prof.'});
htl.FontSize=7;
%set(htl,'fontsize',12)
xlabel('Longitude')
ylabel('Latitude')
title(['Float WMO ' floatnum])
if PARAM.PRINT==1
    figure(hf)
    set(hf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]);
    %eval(['print -depsc2 ' plotpath floatnum '_pos_' titflag titsavedata '.eps']);
    eval(['print -f' num2str(hf.Number) ' -dpng ' plotpath floatnum '_pos_' titflag titsavedata '.png']);
end



% carte 2 (positions)
% -------------------

zone_visu=[floor(min(lat))-10 ceil(max(lat))+10  floor(min(lon))-20 ceil(max(lon))+20];

reso='LR';proj='miller';
[hf,ha]=fct_pltmap(zone_visu,reso,proj);

%h=m_plot(lon,lat,'k-')
[elev,lonbath,latbath]=m_tbase([floor(min(lon))-25 ceil(max(lon))+25 floor(min(lat))-10 ceil(max(lat))+10]);
m_contour(lonbath,latbath,elev,[-4000:1000:0],'color',[0.5 0.5 0.5])
for ii=2:length(lon)-1
    m_plot(lon(ii),lat(ii),'color',map(ii,:),'marker','o','markerfacecolor',map(ii,:),'markersize',8)
end

ii=1;
hl(1)=m_plot(lon(ii),lat(ii),'color',map(ii,:),'marker','d','markerfacecolor',map(ii,:),'markersize',8);
ii=length(lon);
hl(2)=m_plot(lon(ii),lat(ii),'color',map(ii,:),'marker','s','markerfacecolor',map(ii,:),'markersize',8);
%htl=legend(hl,'First prof.','Last prof.','bestinside');
htl=legend(hl,{'First prof.','Last prof.'});

set(htl,'fontsize',12)
xlabel('Longitude')
ylabel('Latitude')
title(['Float WMO ' floatnum])

if PARAM.PRINT==1
    set(gcf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]);
    %eval(['print -depsc2 ' plotpath floatnum '_pos2_' titflag titsavedata '.eps']);
    eval(['print -f' num2str(hf.Number) ' -dpng ' plotpath floatnum '_pos2_' titflag titsavedata '.png']);
end



figure;
orient landscape
set(gca,'fontsize',18)
hold on
for ii=1:length(lon)
    plot(tabtime(ii),cycnum(ii),'color',map(ii,:),'marker','*')
end

vtick=get(gca,'xtick');
vticklabel=datestr(datenum(1950,1,1)+vtick+juld(1),12);
set(gca,'xtick',vtick,'xticklabel',vticklabel);
set(gca,'ylim',[-2 cycnum(end)]);
grid
ylabel('Cycle number')
xlabel('Date')
title(['Float WMO ' floatnum])
if PARAM.PRINT==1
    set(gcf,'papertype','usletter','paperunits','inches','paperorientation','landscape','paperposition',[.25,.75,9.5,8]);
    %eval(['print -depsc2 ' plotpath floatnum '_cycle_' titflag titsavedata '.eps']);
    eval(['print -dpng ' plotpath floatnum '_cycle_' titflag titsavedata '.png']);
end


if isplot_ice>1
    FLm=format_flags_char2num(FLm);
    figure(2)
    %subplot(3,1,3)
    nexttile
    full_cycle=[1:max(FLm.cycle_number.data)];
    empty_cycle=~ismember(full_cycle,FLm.cycle_number.data);
    b=bar(full_cycle,empty_cycle*10);
    b.FaceColor=[0.9 0.9 0.9];
    b.LineStyle='none'
    hold on
    grid on
    box on
    title ('Position Flag QC')
    xlabel('Cycle Number')
    ylabel('Flag')
    %scatter(FLm.longitude.data,FLm.latitude.data,3,FLm.position_qc.data)
    m=max(FLm.position_qc.data);
    plot(FLm.cycle_number.data,FLm.position_qc.data,'.-')
    set(gca,'YLim',[0 m])
    
    if isplot==1
        fileout=[CONFIG.DIR_PLOT 'preliminaire/' floatname '/' floatname '_ice_param.pdf'];
        if exist([CONFIG.DIR_PLOT 'preliminaire/' floatname '/'])==0
            mkdir([CONFIG.DIR_PLOT 'preliminaire/' floatname '/'])
        end
        set(gcf,'Position',[1 46 660 800],'paperPositionMode','auto');
        print(fileout,'-dpdf')
    end
end

%clear all
%close all
disp('plotdata termine');
% Next line commented by T. Reynaud 11/09/2020
%eval('cd ..');

