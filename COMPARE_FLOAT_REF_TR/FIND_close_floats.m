% -========================================================
%   USAGE :   FIND_close_floats(floatname,dacname,varargin)
%   PURPOSE : find all the profiles available on GDAC, close to the profiles of a given float, at a given distance and in a given period of time
%             and compare salinity on Theta levels for each float profile. Anomalies between float profiles and GDAC profiles are computed below PRES_MIN.
% -----------------------------------
%   INPUT :
%    floatname  (char)  e.g. '690258'
%    dacname    (char) e.g.  'coriolis'
%
%   OPTIONNAL INPUT :
%    'ECARLON' (float)  maximum difference in longitude  (default is 0.2 deg)
%    'ECARLAT' (float)  maximum difference in latitude  (default is 0.1 deg)
%    'ECARDAY' (float)   maximum difference in time (default is 365 days)
%    'ZOOM_DEPTH' (float) used to zoom on the deepest level on the plots (default is 0 i.e no zoom)
%    'PRES_MIN' (float)  minimum  float pressure used to compute salinity anomalies (default is 1000db)
%    'EXCLUDE_FLOATS (cell of char)  e.g. {'6900458','6900459'} list of
%    floats to be excluded from the comparison, (default, empty)
%    'OCEAN' (cell of char) 'A' 'M' 'P' or 'I'  ex: {'A','P'} , 'M' is for Mediterranee, 'A' is Atlantic (without Med)
%    'UPDATE'    (logical)      UPDATE=1 if the list of all Argo profiles availble
%                               on GDAC is updated each time this program is run, UPDATE =0 otherwise (default=1)
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
function FIND_close_floats(floatname,dacname,varargin)
close all

CONFIG=load_configuration('config.txt');

% To save plots
CONFIG.plotpathini=[CONFIG.DIR_PLOT '/ref_database/'];
CONFIG.plotpath=[CONFIG.plotpathini '/' floatname '/'];
if exist(CONFIG.plotpath)~=7
    mkdir(CONFIG.plotpath)
end

% To save list of nearby floats
CONFIG.file_nearby_floats=[CONFIG.plotpath '/list_' floatname '.mat'];

n=length(varargin)

if n/2~=floor(n/2)
    error('check the imput arguments')
end

f=varargin(1:2:end)
c=varargin(2:2:end)
s = cell2struct(c,f,2);

% default CONFIG
PARAM.ECARLON=0.2;
PARAM.ECARLAT=0.1;
PARAM.ECARDAY=365;
PARAM.ZOOM_DEPTH=0;
PARAM.PRES_MIN=1000;
PARAM.EXCLUDE_FLOATS={''};
PARAM.OCEAN={'A','M','P','I'};
PARAM.UPDATE=1;
PARAM.PRINT=0;
PARAM.DATATYPE='raw';

% Input CONFIG
if isfield(s,'ECARLON')==1;PARAM.ECARLON=s.ECARLON;end;
if isfield(s,'ECARLAT')==1;PARAM.ECARLAT=s.ECARLAT;end;
if isfield(s,'ECARDAY')==1;PARAM.ECARDAY=s.ECARDAY;end;
if isfield(s,'ZOOM_DEPTH')==1;PARAM.ZOOM_DEPTH=s.ZOOM_DEPTH;end;
if isfield(s,'PRES_MIN')==1;PARAM.PRES_MIN=s.PRES_MIN;end;
if isfield(s,'EXCLUDE_FLOATS')==1;PARAM.EXCLUDE_FLOATS=s.EXCLUDE_FLOATS;end;
if isfield(s,'OCEAN')==1;PARAM.OCEAN=s.OCEAN;end;
if isfield(s,'UPDATE')==1;PARAM.UPDATE=s.UPDATE;end;
if isfield(s,'PRINT')==1;PARAM.PRINT=s.PRINT;end;
if isfield(s,'DATATYPE')==1;PARAM.DATATYPE=s.DATATYPE;end;

NcVar.latitude.name='LATITUDE';
NcVar.longitude.name='LONGITUDE';
F=read_netcdf_allthefile([CONFIG.DIR_FTP dacname '/' floatname '/' floatname '_prof.nc'],NcVar);
time_slot=[1:length(F.latitude.data)];
BOITE.lonmin=floor(min(F.longitude.data(time_slot))*10)/10;
BOITE.lonmax=ceil(max(F.longitude.data(time_slot))*10)/10;
BOITE.latmin=floor(min(F.latitude.data(time_slot))*10)/10;
BOITE.latmax=ceil(max(F.latitude.data(time_slot))*10)/10;
BOITE.name='ATLN';
BOITE.shiftEW='grwch';

plotpathini= CONFIG.plotpathini;

listfloat={floatname}


% flotteurs a exclure de la comparaison
exclude_float=PARAM.EXCLUDE_FLOATS;

% index_file_name=[CONFIG.DIR_FTP '../ar_index_global_prof.txt'];
% fr=fopen(index_file_name,'r')
% index_data = textscan(fr,'%s%s%s%s%s%s%s%s' , 'delimiter' , ',' ,'commentStyle','#','Headerlines',9);
% % transforme les longitue latitude   en double %  (il vaut mieux tout lire en %s, cela evite que textscan plante si ce qui est dans le fichier index ne correspond pas a une valeur numerique)
% index_data{4} = str2double(index_data{4});
% index_data{3} = str2double(index_data{3});
% fclose(fr);
% [liste_profile_files,liste_floats_files] = select_floats_from_index(index_data,'lonmin',BOITE.lonmin, 'lonmax', BOITE.lonmax, 'latmin',BOITE.latmin, 'latmax', BOITE.latmax ,'ocean','A','exclude_type',{'852','851'},'datemin','20020101');

% trouve les flotteurs TR dans la zone
if ~exist(CONFIG.file_nearby_floats) || PARAM.UPDATE==1
    disp('reading GDAC index file ...')
    index_file_name=[CONFIG.DIR_REFDATA_CORIOLIS 'ar_index_global_prof.txt'];
    if ~exist(index_file_name)
        index_file_name=[CONFIG.DIR_REFDATA_CORIOLIS '../ar_index_global_prof.txt'];
    end
    index_data = read_index_prof(index_file_name);
    
    [liste_profile_files,liste_floats_files] = select_floats_from_index(index_data,'lonmin',BOITE.lonmin, 'lonmax', BOITE.lonmax, 'latmin',BOITE.latmin, 'latmax', BOITE.latmax ,'ocean',PARAM.OCEAN,'exclude_type',{'852','851'},'datemin','20020101');
    save(CONFIG.file_nearby_floats,'liste_profile_files','liste_floats_files')
else
    load(CONFIG.file_nearby_floats);
    disp(['Number of  nearby floats selected : ' num2str(length(liste_floats_files))])
end

[tabdac,r]=strtok(liste_floats_files,'/');
tabfloat=strtok(r,'/');


iio=ismember(tabfloat,exclude_float);
tabfloat=tabfloat(~iio);
tabdac=tabdac(~iio);


% addpath('/home1/homedir5/perso/ccabanes/dvlpRD/Argo/Lib_forCathy/ObsInSitu/Lib_Argo/')
% addpath('/home1/homedir5/perso/ccabanes/dvlpRD/Argo/Lib_forCathy/ObsInSitu/Lib_Argo/Plots/')
% addpath(DIR_OW_ORIG) % codes OW originaux
%
% % sauve le path
% path_save=path;

% Topo_ficin='/home5/pharos/argo/DMARGO/OW/TOPO/topo.onetenthdeg.nc';
% Topo=read_netcdf_allthefile(Topo_ficin);

[Topo,DimT]=read_netcdf_allthefile(CONFIG.FILE_TOPO);


NcVar.wmo_inst_type.name='WMO_INST_TYPE';
NcVar.pi_name.name='PI_NAME';
NcVar.float_serial_no.name='FLOAT_SERIAL_NO';
NcVar.latitude.name='LATITUDE';
NcVar.longitude.name='LONGITUDE';
NcVar.juld.name='JULD';
NcVar.cycle_number.name='CYCLE_NUMBER';
NcVar.paltform_number.name='PLATFORM_NUMBER';


iboite_prec=1;
ALL=struct();

%  Lit la liste de tous les flotteurs et stocke les informations necessaires (longitude, latitude, date, cycle)
for ifloat=1:length(tabfloat)
    
    thecal='';
    flt_name=tabfloat{ifloat};
    flt_dir='';
    flt_dir = deblank(flt_dir);
    flt_name = deblank(num2str(flt_name));
    dacname = deblank(tabdac{ifloat});
    
    filename = [CONFIG.DIR_REFDATA_CORIOLIS  strtrim(dacname) '/' strtrim(flt_name) '/' strtrim(flt_name) '_prof.nc'];
    F = read_netcdf_allthefile(filename,NcVar);
    F = shiftEW(F,'longitude',BOITE.shiftEW);
    wmo=str2num(F.wmo_inst_type.data(1,:));
    
    % vérifie si le flotteur passe dans la boite
    iloc=(F.latitude.data<=BOITE.latmax&F.latitude.data>=BOITE.latmin&F.longitude.data>=BOITE.lonmin&F.longitude.data<=BOITE.lonmax);
    if sum(iloc)>0
        length_float=length(F.latitude.data);
        iboite = iboite_prec+length_float-1;
        
        ALL.latitude.data(iboite_prec:iboite)=F.latitude.data;
        ALL.longitude.data(iboite_prec:iboite)=F.longitude.data;
        ALL.juld.data(iboite_prec:iboite)=F.juld.data;
        ALL.cycle_number.data(iboite_prec:iboite)=F.cycle_number.data;
        ALL.floatname.data(iboite_prec:iboite,:)=str2num(char(F.platform_number.data));
        iboite_prec = iboite+1;
    end
end

thefloatlist=str2num(char(listfloat));
iin=~ismember(ALL.floatname.data,thefloatlist);
iou=ismember(ALL.floatname.data,thefloatlist);
thefield=fieldnames(ALL)
for k=1:length(thefield)
    REF.(thefield{k}).data = ALL.(thefield{k}).data(iin);
    LIS.(thefield{k}).data = ALL.(thefield{k}).data(iou);
end
ipr=0;
clear thefloat_tocheck thecycle_tocheck thefloat_neighb thecycle_neighb


for k=1:length(listfloat)
    iou = ismember(LIS.floatname.data,thefloatlist(k));
    latitude = LIS.latitude.data(iou);
    longitude = LIS.longitude.data(iou);
    juld = LIS.juld.data(iou);
    cycle = LIS.cycle_number.data(iou);
    for ij=1:sum(iou)
        iproche = find(abs(REF.latitude.data-latitude(ij))< PARAM.ECARLAT & abs(REF.longitude.data-longitude(ij))< PARAM.ECARLON & abs(REF.juld.data-juld(ij))< PARAM.ECARDAY);
        if length(iproche)>0
            %keyboard
            for ijk=1:length(iproche)
                ipr = ipr+1;
                thefloat_tocheck(ipr) = thefloatlist(k);
                thecycle_tocheck(ipr) = cycle(ij);
                thefloat_neighb(ipr) = REF.floatname.data(iproche(ijk));
                thecycle_neighb(ipr) = REF.cycle_number.data(iproche(ijk));
            end
        end
    end
end


% trace les comparaisons pour chaque flotteur (carte, profils et diag theta/S)
uniq_float_to_check = unique(thefloat_tocheck);

themindepth=PARAM.ZOOM_DEPTH;
clear themean thefloatused
for k=1:length(uniq_float_to_check)
    
    flt_name = num2str(uniq_float_to_check(k));
    plotpath=[plotpathini '/' flt_name '/'];
    if exist(plotpath)~=7
        mkdir(plotpath)
    end
    iidac = ismember(tabfloat,flt_name);
    
    dacname = deblank(tabdac{iidac});
    if strcmp(PARAM.DATATYPE,'adj')
		vertical_sampling_scheme='Primary sampling';
		IncludeDescProf=1;
		[file_list] = select_float_files_on_ftp(floatname,dacname,CONFIG.DIR_FTP,'C',IncludeDescProf);
		[F,Dim,thelist_ext2]=create_multi_from_filelist(floatname,dacname,CONFIG.DIR_FTP,file_list,vertical_sampling_scheme,'');
	else
		 filename=[CONFIG.DIR_FTP  dacname '/' floatname '/' floatname '_prof.nc'];
		 F=read_netcdf_allthefile(filename);
	end
    %filename = [CONFIG.DIR_FTP  strtrim(dacname) '/' strtrim(flt_name) '/' strtrim(flt_name) '_prof.nc'];
    %F = read_netcdf_allthefile(filename);
    F = replace_fill_bynan(F);
    F = format_flags_char2num(F);
	F.tpot=F.temp;
	if strcmp(PARAM.DATATYPE,'adj')
        F.psal=F.psal_adjusted;F.temp=F.temp_adjusted;F.pres=F.pres_adjusted;
		F.psal_qc=F.psal_adjusted_qc;F.temp_qc=F.temp_adjusted_qc;F.pres_qc=F.pres_adjusted_qc;

    end
    F.tpot.data = sw_ptmp(F.psal.data,F.temp.data,F.pres.data,0);
    F=shiftEW(F,'longitude',BOITE.shiftEW);
    
    region.lonmin=min(F.longitude.data)-10;
    region.lonmax=max(F.longitude.data)+10;
    region.latmin=min(F.latitude.data)-5;
    region.latmax=max(F.latitude.data)+5;
    region.shiftEW=BOITE.shiftEW;
    
    map=jet(length(F.cycle_number.data));
    ik = (ismember(thefloat_tocheck,uniq_float_to_check(k)));
    
    unique_cycle_to_check=unique(thecycle_tocheck(ik));
    
    for icycle=1:length(unique_cycle_to_check)
        icycle
        
        theprof=find(F.cycle_number.data==unique_cycle_to_check(icycle)&F.direction.data=='A');
        if isempty(theprof)==0
            figure(icycle)  % une figure pour chaque cycle a comparer
            subplot(3,2,1) % carte
            [thetitle] = plot_traj_glob_sans_bathy(F,F.juld.data);
            hold on
             plot_bathy_glob_nofill(region,Topo)
            %          subplot(2,2,2) % theta S a comparer
            %          hold on
            %          theprof=find(F.cycle_number.data==unique_cycle_to_check(icycle)&F.direction.data=='A');
            %          [thetitle]=plot_theta_S_diag(F,theprof,'y');
            %          subplot(2,2,3) % theta S a comparer
            %          hold on
            %          [thetitle]=plot_profile(F,'tpot','pres',theprof,'Mindepth', themindepth)
            %          subplot(2,2,4) % theta S a comparer
            %          hold on
            %          [thetitlei] = plot_profile(F,'psal','pres',theprof,'Mindepth', themindepth)
            ikcy= ik&(ismember(thecycle_tocheck,unique_cycle_to_check(icycle)));
            % trace les données de ref
            uniq_float_neighb = unique(thefloat_neighb(ikcy));
            jkifin=0
            fulltitle=' ';
            
            for kn=1:length(uniq_float_neighb)
                jkideb=jkifin+1;
                ikn = find(ikcy&ismember(thefloat_neighb,uniq_float_neighb(kn)));
                if length(ikn)>0
                    flt_name = num2str(uniq_float_neighb(kn));
                    %                 if strcmp(flt_name,'3901859')
                    %                     keyboard
                    %                 end
                    iidac = ismember(tabfloat,flt_name);
                    dacname = deblank(tabdac{iidac});
                    filename = [CONFIG.DIR_REFDATA_CORIOLIS  strtrim(dacname) '/' strtrim(flt_name) '/' strtrim(flt_name) '_prof.nc'];
                    FL = read_netcdf_allthefile(filename);
                    FL = replace_fill_bynan(FL);
                    FL = format_flags_char2num(FL);
                    FL=construct_best_param(FL,{'temp','pres','psal'},FL);
                    FL.psal=FL.psal_best;FL.temp=FL.temp_best;FL.pres=FL.pres_best;
                    FL.psal.data(FL.psal_qc.data>1&FL.pres.data<2000)=NaN;
					FL.psal.data(FL.psal_qc.data>3&FL.pres.data>2000)=NaN;
                    FL.tpot=FL.temp;
                    FL.tpot.data = sw_ptmp(FL.psal.data,FL.temp.data,FL.pres.data,0);
                    
                    clear theprofn
                    for iprof=1:length(ikn)
                        theres=find(FL.cycle_number.data==thecycle_neighb(ikn(iprof))&FL.direction.data=='A');
                        if isempty(theres)
                            theres=find(FL.cycle_number.data==thecycle_neighb(ikn(iprof))&FL.direction.data=='D');
                        end
                        theprofn(iprof)=theres;
                    end
                    theprofn=unique(theprofn);
                    %                  if length(theprofn)==1
                    %                  theprofn=[theprofn-3:theprofn+3];
                    %                  theprofn(theprofn<1|theprofn>length(FL.latitude.data))=[];
                    %                  end
                    FL=shiftEW(FL,'longitude',BOITE.shiftEW);
                    
                    subplot(3,2,1) % carte
                    hold on
                    %plot_bathy_glob_nofill(region)
                    scatter(FL.longitude.data(theprofn),FL.latitude.data(theprofn),30,'md','filled')
                    % plot_traj_glob_sans_bathy(FL,'bd');
                    
                    subplot(2,2,2) % theta S a comparer
                    hold on
                    
                    [thetitle]=plot_theta_S_diag(FL,theprofn,'n','LineColor','m');
                    
                    subplot(2,2,3) % theta S a comparer
                    hold on
                    [thetitle]=plot_profile(FL,'temp','pres',theprofn,'LineColor','m','Mindepth', themindepth)
                    
                    
                    subplot(2,2,4) % theta S a comparer
                    hold on
                    
                    % [thetitle]=plot_profile(FL,'psal','pres',theprofn,'LineColor','m','Mindepth', themindepth)
                    % keyboard
                    % Difference de salinité sur les niveau theta
                    jkifin=jkideb+length(theprofn)-1
                    for jki=1:length(theprofn)
                        jkicompt=jkideb+jki-1
                        [S_h,P_h]=interp_climatology(FL.psal.data(theprofn(jki),:)',FL.tpot.data(theprofn(jki),:)',FL.pres.data(theprofn(jki),:)',F.psal.data(theprof,:),F.tpot.data(theprof,:),F.pres.data(theprof,:)); % routine OW
                        difference_psal_theta=[F.psal.data(theprof,:)-S_h'];
                        difference_pres_theta=[F.pres.data(theprof,:)-P_h'];
                        ip=find(abs(difference_pres_theta)<150);
                        ipn=find(abs(difference_pres_theta)>150);
                        difference_psal_theta(ipn)=NaN;
                        
                        kk=find(F.pres.data(theprof,ip)>PARAM.PRES_MIN );
                        themean(unique_cycle_to_check(icycle),jkicompt)=mean(difference_psal_theta(ip(kk)));
                        thefloatused(unique_cycle_to_check(icycle),jkicompt)=uniq_float_neighb(kn);
                    end
                    %[thetitle]=plot_profile(FL,'psal','pres',theprofn,'LineColor','m','Mindepth', themindepth)
                    sav=F.psal.data(theprof,:);
                    F.delta_psal=F.psal;
                    F.delta_psal.data(theprof,:)= F.psal.data(theprof,:)-S_h';
                    [thetitlei] = plot_profile(F,'delta_psal','pres',theprof,'Mindepth', themindepth)
                    F.psal.data(theprof,:)=sav;
                    if length(uniq_float_neighb)>1&kn~=length(uniq_float_neighb)
                        fulltitle=[fulltitle thetitle ' ; '];
                    else
                        fulltitle=[fulltitle thetitle];
                    end
                    %title(['Compared to ' thetitle])
                    %keyboard
                end
            end
            
            
            subplot(2,2,2) % theta S a comparer
            %theprof=find(F.cycle_number.data==unique_cycle_to_check(icycle)&F.direction.data=='A');
            [thetitle]=plot_theta_S_diag(F,theprof,'y');
            subplot(2,2,3) % theta S a comparer
            [thetitle]=plot_profile(F,'temp','pres',theprof,'Mindepth', themindepth)
            r=fulltitle;
            ititle=1;
            finaltitle{1} ='Compared to :';
            while isempty(r)==0
                ititle=ititle+1;
                [finaltitle{ititle},r]=strtok(r,';');
            end
            title(finaltitle,'FontSize',8)
            subplot(2,2,4) % theta S a comparer
            hold on
            % F.psal.data(theprof,:)= F.psal.data(theprof,:)-S_h';
            %[thetitlei] = plot_profile(F,'psal','pres',theprof,'Mindepth', themindepth)
            %keyboard
            themean(themean==0)=NaN;
            mean_cycle=meanoutnan(themean(unique_cycle_to_check(icycle)));
            thetitlei=[thetitlei '. Difference PSAL <' num2str(PARAM.PRES_MIN) 'm = ' num2str(mean_cycle)];
            suptitle(thetitlei)
            
            titplt=[num2str(uniq_float_to_check(k)) '_cycle_' num2str(unique_cycle_to_check(icycle)) '.png'];
            eval(['print -dpng ' plotpath  titplt]);
        end
    end
    
    themean(themean==0)=NaN;
    figure(length(unique_cycle_to_check)+1)
    plot(meanoutnan(themean,2),'+')
    hold on
    a=repmat([1:size(themean,1)]',1,size(themean,2));
    b=unique(thefloatused(~thefloatused==0));
    thefloatused_ind=0*thefloatused;
    for llo=1:length(b)
        thefloatused_ind(thefloatused==b(llo))=llo;
    end
    for ijk=1:size(a,2)
        scatter(a(:,ijk),themean(:,ijk),30,thefloatused_ind(:,ijk))
        hold on
    end
    plot(meanoutnan(themean,2),'+')
    grid on
    box on
    title([num2str(uniq_float_to_check(k)) ': comparison with closest RT Argo profiles <' num2str(PARAM.PRES_MIN) 'm'])
    xlabel('cycles')
    ylabel('salinity differences')
    titplt=[num2str(uniq_float_to_check(k)) '_themean1500.png'];
    eval(['print -dpng ' plotpath  titplt]);
    % close all
end
